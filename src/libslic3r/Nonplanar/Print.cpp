#include "Print.hpp"
#include "../I18N.hpp"
#include "../format.hpp"
#include "../Print.hpp"
#include "../Model.hpp"
#include "../PrintConfig.hpp"
#include "../GCodeWriter.hpp"
#include "../GCode.hpp"
#include "../GCode/GCodeProcessor.hpp"
#include "../PlaceholderParser.hpp"
#include "../TriangleMesh.hpp"
#include "../BoundingBox.hpp"
#include "../Polygon.hpp"
#include "../Point.hpp"
#include "../Exception.hpp"
#include <Eigen/Geometry>
#include <boost/nowide/fstream.hpp>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace Slic3r::Nonplanar {
namespace {
[[noreturn]] void invalid(const std::string &s) { throw std::invalid_argument(Slic3r::format(_u8L("Nonplanar slicing: %1%"),s)); }
constexpr const char *placement_key = "orca_nonplanar_placement_v1";

// Stable across project save/reload (object IDs are intentionally excluded).
std::string placement(const Model &model)
{
    nlohmann::json objects = nlohmann::json::array();
    for (const auto *o : model.objects) {
        nlohmann::json object;
        object["printable"] = o->printable;
        for (const auto *i : o->instances) {
            nlohmann::json m = nlohmann::json::array();
            // 3MF writes transforms with finite decimal precision.
            for (int n = 0; n < 16; ++n) m.push_back(std::round(i->get_matrix().matrix().data()[n]*1e5)/1e5);
            object["instances"].push_back(m);
        }
        for (const auto *v : o->volumes) {
            nlohmann::json m = nlohmann::json::array();
            for (int n = 0; n < 16; ++n) m.push_back(std::round(v->get_matrix().matrix().data()[n]*1e5)/1e5);
            object["volumes"].push_back({m,v->mesh().its.vertices.size(),v->mesh().its.indices.size()});
        }
        objects.push_back(object);
    }
    return objects.dump();
}

double scalar(const DynamicPrintConfig &c, const char *key, double fallback)
{
    const auto *o = c.option(key);
    if (!o) return fallback;
    if (const auto *v = dynamic_cast<const ConfigOptionFloats*>(o)) return v->values.front();
    if (const auto *v = dynamic_cast<const ConfigOptionFloatsNullable*>(o)) {
        const double value = v->values.front();
        return std::isfinite(value) ? value : fallback;
    }
    if (const auto *v = dynamic_cast<const ConfigOptionInts*>(o)) return v->values.front();
    return o->type() == coInt ? o->getInt() : o->getFloat();
}
double first(const DynamicPrintConfig &c, const char *key, double fallback)
{
    const auto *o = c.option<ConfigOptionFloats>(key);
    return o && !o->values.empty() ? o->values.front() : fallback;
}
const char *role_name(Role role)
{
    switch (role) {
    case Role::OuterWall: return "Outer wall";
    case Role::InnerWall: return "Inner wall";
    case Role::SparseInfill: return "Sparse infill";
    case Role::TopSkin: return "Top surface";
    case Role::BottomSkin: return "Bottom surface";
    case Role::SolidInfill: return "Internal solid infill";
    }
    return "Internal solid infill";
}
Vec3d vec(XYZ p) { return {p.x,p.y,p.z}; }
} // namespace

std::string source(const Model &model, int plate)
{
    if (!model.model_info) return {};
    const auto &items = model.model_info->metadata_items;
    const auto it = items.find(metadata_key);
    if (it == items.end()) return {};
    const auto version = items.find(version_key);
    if (version == items.end() || version->second != "1") invalid(_u8L("unsupported or missing embedded nonplanar schema version"));
    const auto stack = LayerStack::parse(it->second);
    return stack.plate_index == plate ? stack.source : std::string();
}

void attach(Model &model, const LayerStack &stack)
{
    model.model_info = model.model_info ? std::make_shared<ModelInfo>(*model.model_info) : std::make_shared<ModelInfo>();
    auto &items = model.model_info->metadata_items;
    items[metadata_key] = stack.source;
    items[version_key] = "1";
    items[placement_key] = placement(model);
}
void detach(Model &model)
{
    if (!model.model_info) return;
    model.model_info = std::make_shared<ModelInfo>(*model.model_info);
    for (auto key : {metadata_key,version_key,placement_key}) model.model_info->metadata_items.erase(key);
}

std::shared_ptr<Job> prepare(Print &print, const std::string &text)
{
    auto job = std::make_shared<Job>();
    job->stack = LayerStack::parse(text);
    // Match achievable G-code precision without rejecting finer requests.
    job->stack.limits.geometry_tolerance = std::max(.0002,job->stack.limits.geometry_tolerance);
    job->stack.limits.chord_tolerance = std::max(.0002,job->stack.limits.chord_tolerance);
    job->stack.limits.max_segment_length = std::max(.0004,job->stack.limits.max_segment_length);
    const auto &config = print.full_print_config();
    Settings &s = job->settings;
    s.nozzle_diameter = first(config,"nozzle_diameter",.4);
    const double filament_diameter = first(config,"filament_diameter",1.75), flow_ratio = first(config,"filament_flow_ratio",1.);
    if (!std::isfinite(filament_diameter) || !std::isfinite(flow_ratio) || filament_diameter <= 0. || flow_ratio <= 0.)
        invalid(_u8L("Filament diameter and flow ratio must be positive and finite"));
    s.width = config.get_abs_value("line_width",s.nozzle_diameter);
    if (s.width == 0.) s.width = s.nozzle_diameter*1.125;
    s.walls = config.opt_int("wall_loops");
    s.top_layers = config.opt_int("top_shell_layers");
    s.bottom_layers = config.opt_int("bottom_shell_layers");
    s.density = scalar(config,"sparse_infill_density",15.)*.01;
    s.infill_angle = scalar(config,"infill_direction",45.);
    s.max_z = scalar(config,"printable_height",250.);
    std::vector<Object> objects;
    const Vec3d origin = print.get_plate_origin();
    for (const auto *po : print.objects()) {
        const auto *mo = po->model_object();
        if (mo->volumes.size() != 1) invalid(_u8L("use one watertight model volume per object"));
        for (const auto &instance : po->instances()) {
            Object object;
            object.id = po->id().id;
            for (const auto *volume : mo->volumes) {
                if (!volume->is_model_part()) invalid(_u8L("modifier and negative volumes are unsupported; use a watertight mesh with native cavities"));
                const Transform3d transform = instance.model_instance->get_matrix()*volume->get_matrix();
                const auto &mesh = volume->mesh().its;
                for (const auto &face : mesh.indices) {
                    Triangle triangle;
                    for (int n = 0; n < 3; ++n) {
                        const Vec3d p = transform*mesh.vertices[face[n]].cast<double>()-origin;
                        triangle.vertices[n] = {p.x(),p.y(),p.z()};
                    }
                    if (transform.linear().determinant() < 0) std::swap(triangle.vertices[1],triangle.vertices[2]);
                    object.triangles.push_back(triangle);
                }
            }
            objects.push_back(std::move(object));
        }
    }
    job->layers = slice(job->stack,objects,s,[&print] { if (print.canceled()) throw CanceledException(); });
    return job;
}

std::string export_gcode(Print &print, const Job &job, const std::string &path, GCodeProcessorResult *result)
{
    const auto &config = print.full_print_config();
    GCodeWriter writer;
    writer.apply_print_config(print.config());
    writer.config.high_precision_gcode.value = true;
    writer.set_is_bbl_machine(print.is_BBL_printer());
    writer.set_extruders({0});
    writer.set_extruder(0);
    auto parser = print.placeholder_parser();
    parser.set("initial_tool",0);
    parser.set("initial_extruder",0);
    for (const char *key : {"initial_no_support_tool","initial_no_support_extruder","initial_no_support_hotend",
                           "initial_filament_id","initial_no_support_filament_id","initial_nozzle_id",
                           "current_filament_id","current_extruder_id","current_nozzle_id","current_hotend"})
        parser.set(key,0);
    parser.set("current_extruder",0);
    // Printer scripts still need physical-tool vectors, including unused tools.
    std::vector<int> first_filaments(print.config().nozzle_diameter.size(),-1);
    first_filaments.front() = 0;
    parser.set("first_tools",new ConfigOptionInts(first_filaments));
    parser.set("first_filaments",new ConfigOptionInts(first_filaments));
    parser.set("first_non_support_filaments",new ConfigOptionInts(first_filaments));
    parser.set("nozzle_diameter_at_nozzle_id",new ConfigOptionFloats{job.settings.nozzle_diameter});
    parser.set("overall_chamber_temperature",config.option<ConfigOptionInts>("chamber_temperature")->values.front());
    parser.set("total_layer_count",int(job.layers.size()));
    size_t last_printed = 0, first_printed = 0;
    const Path *finishing_path = nullptr;
    double max_printed_z = 0.;
    for (size_t k = 0; k < job.layers.size(); ++k)
        for (const Path &path : job.layers[k].paths)
            if (std::any_of(path.segments.begin(),path.segments.end(),[](const Segment &s) { return s.volume > 0.; })) {
                if (!finishing_path) first_printed = k;
                last_printed = k; finishing_path = &path;
                max_printed_z = std::max(max_printed_z,job.layers[k].max_z);
            }
    parser.set("layer_num",0);
    parser.set("layer_z",job.layers[first_printed].max_z);
    parser.set("max_layer_z",max_printed_z);
    parser.set("max_print_z",std::ceil(max_printed_z));
    double first_layer_area = 0.;
    for (const auto &path : job.layers[first_printed].paths) first_layer_area += path.projected_area;
    parser.set("hold_chamber_temp_for_flat_print",max_printed_z > 0. && max_printed_z < .3 && first_layer_area > 40000.);
    parser.set("max_print_height",job.settings.max_z);
    parser.set("z_offset",0.);
    parser.set("model_name",print.get_model_name());
    parser.set("plate_number",print.get_plate_number_formatted());
    BoundingBoxf first_layer_box;
    for (const auto &path : job.layers[first_printed].paths) {
        first_layer_box.merge(Vec2d(path.start.x,path.start.y));
        for (const auto &segment : path.segments) first_layer_box.merge(Vec2d(segment.end.x,segment.end.y));
    }
    first_layer_box.offset(job.settings.width*.5);
    const Vec2d size = first_layer_box.size(), center = first_layer_box.center();
    parser.set("first_layer_print_min",new ConfigOptionFloats{first_layer_box.min.x(),first_layer_box.min.y()});
    parser.set("first_layer_print_max",new ConfigOptionFloats{first_layer_box.max.x(),first_layer_box.max.y()});
    parser.set("first_layer_print_size",new ConfigOptionFloats{size.x(),size.y()});
    parser.set("first_layer_center_no_wipe_tower",new ConfigOptionFloats{center.x(),center.y()});
    const BedType bed_type = config.opt_enum<BedType>("curr_bed_type");
    const auto *initial_bed = config.option<ConfigOptionInts>(get_bed_temp_1st_layer_key(bed_type));
    const auto *bed = config.option<ConfigOptionInts>(get_bed_temp_key(bed_type));
    const auto *initial_nozzle = config.option<ConfigOptionInts>("nozzle_temperature_initial_layer");
    const auto *nozzle = config.option<ConfigOptionInts>("nozzle_temperature");
    parser.set("bed_temperature_initial_layer",initial_bed->clone());
    parser.set("first_layer_bed_temperature",initial_bed->clone());
    parser.set("bed_temperature",bed->clone());
    parser.set("first_layer_temperature",initial_nozzle->clone());
    parser.set("bed_temperature_initial_layer_single",initial_bed->values.front());
    double flush_speed = first(config,"filament_flush_volumetric_speed",0.);
    if (flush_speed == 0.) flush_speed = first(config,"filament_max_volumetric_speed",0.);
    const bool fast_flush = config.opt_enum<PrimeVolumeMode>("prime_volume_mode") == PrimeVolumeMode::pvmFast;
    const auto *flush_config = config.option<ConfigOptionInts>(fast_flush ? "filament_flush_temp_fast" : "filament_flush_temp");
    int flush_temperature = flush_config->values.front();
    if (flush_temperature == 0) flush_temperature = config.option<ConfigOptionInts>("nozzle_temperature_range_high")->values.front();
    parser.set("flush_volumetric_speeds",new ConfigOptionFloats{flush_speed});
    parser.set("flush_temperatures",new ConfigOptionInts{flush_temperature});
    const double initial_height = job.stack.vertical_gap(first_printed,XY{job.layers[first_printed].paths.front().start.x,job.layers[first_printed].paths.front().start.y});
    const double flow_limit = first(config,"filament_max_volumetric_speed",0.);
    double outer_flow = scalar(config,"outer_wall_speed",40.)*job.settings.width*initial_height;
    if (flow_limit > 0.) outer_flow = std::min(outer_flow,flow_limit);
    parser.set("outer_wall_volumetric_speed",outer_flow);
    parser.set("initial_layer_print_height",initial_height);
    parser.set("during_print_exhaust_fan_speed_num",config.option("during_print_exhaust_fan_speed")->clone());
    PlaceholderParser::ContextData context;
    const std::string startup = parser.process(config.opt_string("machine_start_gcode"),0,nullptr,&context);
    const auto &filament_start_template = config.option<ConfigOptionStrings>("filament_start_gcode")->values.front();
    const auto &filament_end_template = config.option<ConfigOptionStrings>("filament_end_gcode")->values.front();
    const std::string filament_start = parser.process(filament_start_template,0,nullptr,&context);
    parser.set("layer_num",int(last_printed));
    parser.set("layer_z",finishing_path->segments.back().end.z);
    const std::string shutdown = parser.process(config.opt_string("machine_end_gcode"),0,nullptr,&context);
    const std::string filament_end = parser.process(filament_end_template,0,nullptr,&context);
    const double diameter = first(config,"filament_diameter",1.75);
    const double filament_area = .25*std::acos(-1.)*diameter*diameter;
    const double ratio = first(config,"filament_flow_ratio",1.);
    const double max_flow = first(config,"filament_max_volumetric_speed",0.);
    const double speed = scalar(config,"outer_wall_speed",40.);
    const double travel = scalar(config,"travel_speed",120.);
    std::array<double,4> axis_limits {first(config,"machine_max_speed_x",0.),first(config,"machine_max_speed_y",0.),
                                      first(config,"machine_max_speed_z",0.),first(config,"machine_max_speed_e",0.)};
    double xy_travel = travel;
    for (size_t axis = 0; axis < 2; ++axis) if (axis_limits[axis] > 0.) xy_travel = std::min(xy_travel,axis_limits[axis]);
    const double z_travel = axis_limits[2] > 0. ? std::min(travel,axis_limits[2]) : travel;
    writer.config.travel_speed.values = {xy_travel};
    writer.config.initial_layer_travel_speed.values = {FloatOrPercent(xy_travel,false)};
    writer.config.travel_speed_z.values = {z_travel};
    const double tol = job.stack.limits.geometry_tolerance;
    auto segment_feed = [&](const Path &path, const Segment &segment, const Vec3d &previous, size_t k) {
        const char *key = path.role == Role::OuterWall ? "outer_wall_speed" :
                          path.role == Role::InnerWall ? "inner_wall_speed" :
                          path.role == Role::SparseInfill ? "sparse_infill_speed" :
                          path.role == Role::TopSkin ? "top_surface_speed" : "internal_solid_infill_speed";
        double requested = scalar(config,key,speed);
        if (requested == 0.) requested = speed; // Orca's automatic speed, conservatively bounded.
        if (k == 0) requested = std::min(requested,scalar(config,
            path.role == Role::OuterWall || path.role == Role::InnerWall ? "initial_layer_speed" : "initial_layer_infill_speed",requested));
        if (requested <= 0. || !std::isfinite(requested)) invalid(_u8L("native deposition speed must be positive and finite"));
        const Vec3d end = vec(segment.end);
        const double length = (end-previous).norm(), rate = segment.volume*ratio/length;
        double feed = max_flow > 0 ? std::min(requested,max_flow/rate) : requested;
        for (size_t axis = 0; axis < 3; ++axis) {
            const double fraction = std::abs(end[axis]-previous[axis])/length;
            if (fraction > 0. && axis_limits[axis] > 0.) feed = std::min(feed,axis_limits[axis]/fraction);
        }
        if (axis_limits[3] > 0.) feed = std::min(feed,axis_limits[3]*filament_area/rate);
        return feed;
    };
    // Build in memory before opening the output: validation errors never leave
    // a nominally successful, partly generated machine file.
    std::ostringstream out;
    out.imbue(std::locale::classic()); out << std::fixed << std::setprecision(6);
    out << "; generated by OrcaSlicer native nonplanar\n; orca_nonplanar_schema_version = 1\n"
        << "; orca_nonplanar_layers = " << job.stack.source << "\n"
        << "; PREVIEW_WIDTH: " << job.settings.width << "\n; PREVIEW_HEIGHT: " << scalar(config,"layer_height",.2)
        << "\n; FEATURE: Custom\n";
    out << writer.set_bed_temperature(initial_bed->values.front(),false)
        << writer.set_temperature(initial_nozzle->values.front(),false)
        << startup << '\n'
        << writer.set_bed_temperature(initial_bed->values.front(),true)
        << writer.set_temperature(initial_nozzle->values.front(),true)
        << filament_start << '\n' << writer.preamble() << writer.reset_e(true);
    double printed_max = job.stack.bounds(job.stack.bed).second;
    double volume = 0.;
    int64_t previous_feed = -1;
    auto clear_travel = [&](const Vec3d &target) {
        const double safe_z = std::min(job.settings.max_z,
            std::max(printed_max,target.z())+std::max(1.,job.stack.limits.max_gap)+tol);
        out << writer.retract()
            << writer.travel_to_z(safe_z,"nonplanar clearance",true)
            << writer.travel_to_xyz({target.x(),target.y(),safe_z},"nonplanar travel",true)
            << writer.travel_to_xyz(target,"nonplanar approach",true) << writer.unretract();
        // Travel and retraction emit their own modal feed rates.
        previous_feed = -1;
    };
    for (size_t k = 0; k < job.layers.size(); ++k) {
        if (print.canceled()) throw CanceledException();
        const Layer &layer = job.layers[k];
        if (k == 1) out << writer.set_bed_temperature(bed->values.front(),false) << writer.set_temperature(nozzle->values.front(),false);
        parser.set("layer_num",int(k));
        parser.set("layer_z",layer.max_z);
        out << "; FEATURE: Custom\n";
        for (const char *key : {"before_layer_change_gcode","layer_change_gcode","time_lapse_gcode"})
            if (config.has(key) && !config.opt_string(key).empty())
                out << parser.process(config.opt_string(key),0,nullptr,&context) << '\n';
        out << writer.preamble() << writer.reset_e(true);
        out << "; CHANGE_LAYER\n; Z_HEIGHT: " << layer.max_z << '\n';
        double seconds = 0.;
        for (const auto &path : layer.paths) {
            Vec3d previous = vec(path.start);
            for (const auto &segment : path.segments) {
                if (segment.volume > 0.) seconds += (vec(segment.end)-previous).norm()/segment_feed(path,segment,previous,k);
                previous = vec(segment.end);
            }
        }
        const double cooling_time = first(config,"slow_down_layer_time",0.);
        const bool slow_down = config.option<ConfigOptionBools>("slow_down_for_layer_cooling")->values.front();
        const double cooling_scale = slow_down && cooling_time > 0. ? std::min(1.,seconds/cooling_time) : 1.;
        const double fan_time = first(config,"fan_cooling_layer_time",cooling_time);
        const double fan_fraction = fan_time > cooling_time ? std::clamp((fan_time-seconds)/(fan_time-cooling_time),0.,1.) : 1.;
        double fan = first(config,"fan_min_speed",0.) + fan_fraction*(first(config,"fan_max_speed",0.)-first(config,"fan_min_speed",0.));
        const int off_layers = config.option<ConfigOptionInts>("close_fan_the_first_x_layers")->values.front();
        const int full_layer = config.option<ConfigOptionInts>("full_fan_speed_layer")->values.front();
        if (int(k) < off_layers) fan = 0.;
        else if (k == 0 && config.option<ConfigOptionInts>("initial_layer_fan_speed")->values.front() >= 0)
            fan = config.option<ConfigOptionInts>("initial_layer_fan_speed")->values.front();
        else if (full_layer > off_layers && int(k+1) < full_layer)
            fan *= double(int(k+1)-off_layers)/(full_layer-off_layers);
        out << writer.set_fan(unsigned(std::lround(fan)));
        for (const auto &path_data : layer.paths) {
            if (std::none_of(path_data.segments.begin(),path_data.segments.end(),[](const Segment &s) { return s.volume > 0.; })) continue;
            const Vec3d start = vec(path_data.start);
            // The available toolhead profile is a radius and a bottom clearance
            // height. Conservatively clear the whole printed height envelope.
            clear_travel(start);
            out << "; FEATURE: " << role_name(path_data.role) << '\n';
            Vec3d previous = start;
            for (const Segment &segment : path_data.segments) {
                const Vec3d end = vec(segment.end);
                if (segment.volume == 0.) { clear_travel(end); previous = end; continue; }
                const double maximum_feed = segment_feed(path_data,segment,previous,k);
                const double feed = std::min(maximum_feed,std::max(first(config,"slow_down_min_speed",0.),maximum_feed*cooling_scale));
                std::string move = writer.extrude_to_xyz(end,segment.volume*ratio/filament_area);
                const double F = feed*60.;
                const int64_t formatted_feed = std::llround(F*GCodeFormatter::pow_10[GCodeFormatter::F_EXPORT_DIGITS]);
                if (formatted_feed != previous_feed) {
                    // Reuse the writer's formatting and speed state, but put
                    // the changed F word on this move instead of a separate G1.
                    std::string speed = writer.set_speed(F);
                    speed.erase(0,2); // G1
                    speed.pop_back(); // newline
                    move.insert(2,speed);
                    previous_feed = formatted_feed;
                }
                out << move;
                printed_max = std::max(printed_max,segment.max_top_z); volume += segment.volume*ratio;
                previous = end;
            }
        }
    }
    out << writer.retract() << "; FEATURE: Custom\n" << filament_end << '\n' << shutdown << '\n' << writer.postamble();
    std::string full_config;
    GCode::append_full_config(print,full_config);
    out << "; CONFIG_BLOCK_START\n" << full_config << "; CONFIG_BLOCK_END\n";
    boost::nowide::ofstream file(path,std::ios::binary);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << out.str(); file.close();
    auto &statistics = print.print_statistics();
    statistics.clear(); statistics.total_extruded_volume = volume;
    statistics.total_used_filament = volume/filament_area;
    statistics.filament_stats[0] = statistics.total_used_filament;
    statistics.total_weight = volume*.001*first(config,"filament_density",1.24);
    statistics.total_cost = statistics.total_weight*.001*first(config,"filament_cost",0.);
    if (result) {
        GCodeProcessor processor;
        processor.apply_config(print.config());
        processor.process_file(path,[&print] { if (print.canceled()) throw CanceledException(); });
        statistics.estimated_normal_print_time = processor.get_time_dhm(PrintEstimatedStatistics::ETimeMode::Normal);
        *result = processor.extract_result();
    }
    return path;
}
} // namespace Slic3r::Nonplanar
