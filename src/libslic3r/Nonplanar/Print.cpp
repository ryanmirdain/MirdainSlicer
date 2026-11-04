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
#include <unordered_map>
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
        if (!std::isfinite(value)) invalid(Slic3r::format(_u8L("%1% must be explicitly set for nonplanar slicing"),key));
        return value;
    }
    if (const auto *v = dynamic_cast<const ConfigOptionInts*>(o)) return v->values.front();
    return o->type() == coInt ? o->getInt() : o->getFloat();
}
double first(const DynamicPrintConfig &c, const char *key, double fallback)
{
    const auto *o = c.option<ConfigOptionFloats>(key);
    return o && !o->values.empty() ? o->values.front() : fallback;
}
void reject_bool(const DynamicPrintConfig &c, const char *key)
{
    const auto *option = c.option(key);
    if (!option) return;
    const auto *values = dynamic_cast<const ConfigOptionBools*>(option);
    const bool enabled = values ? std::any_of(values->values.begin(),values->values.end(),[](auto v) { return v != 0; }) : option->getBool();
    if (enabled) invalid(Slic3r::format(_u8L("%1% is unsupported; disable it explicitly"),key));
}
void validate_config(const DynamicPrintConfig &c)
{
    for (const auto *key : {"enable_support","enable_prime_tower","spiral_mode","belt_printer",
                           "enable_overhang_speed", "enable_wrapping_detection", "detect_thin_wall"})
        reject_bool(c,key);
    if (c.opt_enum<ArcFittingMode>("enable_arc_fitting") != ArcFittingMode::Disabled)
        invalid(_u8L("arc fitting is unsupported; disable it explicitly"));
    if (c.opt_serialize("sparse_infill_pattern") != "rectilinear") invalid(_u8L("sparse_infill_pattern must be rectilinear"));
    if (c.opt_serialize("wall_generator") != "classic") invalid(_u8L("wall_generator must be classic"));
    if (c.opt_serialize("seam_position") != "aligned" || c.opt_serialize("seam_slope_type") != "none")
        invalid(_u8L("nonplanar walls require aligned seams without scarf joints"));
    reject_bool(c,"staggered_inner_seams");
    for (auto key : {"xy_hole_compensation","xy_contour_compensation","elefant_foot_compensation"})
        if (scalar(c,key,0.) != 0.) invalid(Slic3r::format(_u8L("%1% must be zero to preserve original-model containment"),key));
    if (first(c,"filament_shrink",100.) != 100. || first(c,"filament_shrinkage_compensation_z",100.) != 100.)
        invalid(_u8L("filament shrink compensation must be 100%; the supplied surfaces use original physical geometry"));
    if (c.opt_serialize("gcode_remap_x") != "pos_x" || c.opt_serialize("gcode_remap_y") != "pos_y" || c.opt_serialize("gcode_remap_z") != "pos_z")
        invalid(_u8L("G-code axis remapping is unsupported in nonplanar mode"));
    for (auto key : {"top_surface_pattern","bottom_surface_pattern","internal_solid_infill_pattern"})
        if (c.opt_serialize(key) != "rectilinear") invalid(Slic3r::format(_u8L("%1% must be rectilinear"),key));
    reject_bool(c,"infill_combination");
    if (scalar(c,"top_surface_density",100.) != 100. || scalar(c,"bottom_surface_density",100.) != 100.)
        invalid(_u8L("nonplanar top and bottom skins require 100% density"));
    if (scalar(c,"top_shell_thickness",0.) != 0. || scalar(c,"bottom_shell_thickness",0.) != 0.)
        invalid(_u8L("minimum shell-thickness overrides are unsupported; set them to zero and select top/bottom shell layer counts"));
    if (c.opt_serialize("ironing_type") != "no ironing") invalid(_u8L("ironing is unsupported"));
    if (c.opt_serialize("fuzzy_skin") != "none") invalid(_u8L("fuzzy skin is unsupported"));
    if (c.opt_serialize("print_sequence") != "by layer") invalid(_u8L("print_sequence must be by layer"));
    if (scalar(c,"brim_width",0.) != 0. || scalar(c,"skirt_loops",0.) != 0. || scalar(c,"raft_layers",0.) != 0.)
        invalid(_u8L("brim, skirt and raft must be disabled"));
    if (c.option<ConfigOptionFloats>("nozzle_diameter")->values.size() != 1)
        invalid(_u8L("only single-tool machines are supported"));
    if (c.option<ConfigOptionFloats>("filament_diameter")->values.size() != 1)
        invalid(_u8L("exactly one filament is supported"));
    if (c.opt_serialize("printer_structure") != "undefine" &&
        c.opt_serialize("printer_structure") != "corexy" && c.opt_serialize("printer_structure") != "i3")
        invalid(_u8L("only Cartesian and CoreXY machines are supported"));
    if (c.has("slicing_pipeline_plugin") && !c.option<ConfigOptionStrings>("slicing_pipeline_plugin")->values.empty())
        invalid(_u8L("slicing pipeline plugins are unsupported"));
    if (c.has("post_process") && !c.option<ConfigOptionStrings>("post_process")->values.empty())
        invalid(_u8L("postprocessing scripts are unsupported"));
    if (has_bed_exclude_volumes(c) || !c.option<ConfigOptionPoints>("bed_exclude_area")->values.empty())
        invalid(_u8L("machine exclusion regions are not yet supported by nonplanar travel planning"));
    if (scalar(c,"z_offset",0.) != 0.) invalid(_u8L("z_offset must be zero; supply heights in physical plate coordinates"));
    for (auto key : {"before_layer_change_gcode","layer_change_gcode","time_lapse_gcode"})
        if (c.has(key) && !c.opt_string(key).empty()) invalid(Slic3r::format(_u8L("%1% is unsupported in nonplanar mode"),key));
}

void validate_mesh(const indexed_triangle_set &mesh)
{
    std::unordered_map<uint64_t,std::pair<unsigned,int>> edges;
    for (const auto &face : mesh.indices) for (int n = 0; n < 3; ++n) {
        const uint32_t a = face[n], b = face[(n+1)%3];
        const uint64_t key = (uint64_t(std::min(a,b)) << 32) | std::max(a,b);
        auto &edge = edges[key]; ++edge.first; edge.second += a < b ? 1 : -1;
    }
    for (const auto &entry : edges)
        if (entry.second.first != 2 || entry.second.second != 0)
            invalid(_u8L("original mesh must be closed, manifold and consistently oriented"));
    if (its_volume(mesh) <= 0.) invalid(_u8L("original mesh must have outward orientation and positive volume"));
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
    if (job->stack.limits.geometry_tolerance < .0002 || job->stack.limits.chord_tolerance < .0002 ||
        job->stack.limits.max_segment_length <= .0002)
        invalid(_u8L("native G-code uses four decimal XYZ coordinates; geometry/chord tolerances must be >= 0.0002 mm and max segment length > 0.0002 mm"));
    const auto &model = print.model();
    const auto &items = model.model_info->metadata_items;
    const auto saved = items.find(placement_key);
    if (saved == items.end() || saved->second != placement(model))
        invalid(_u8L("model placement changed since layer import; restore placement or import redesigned layers"));
    const auto &config = print.full_print_config();
    validate_config(config);
    if (!model.get_curr_plate_custom_gcodes().gcodes.empty()) invalid(_u8L("per-height custom G-code is unsupported"));
    Settings &s = job->settings;
    s.nozzle_diameter = first(config,"nozzle_diameter",.4);
    // A zero minimum in supplied layers opts out of preset thickness bounds.
    if (job->stack.limits.min_gap > 0.) {
        s.min_normal_height = first(config,"min_layer_height",0.);
        s.max_normal_height = first(config,"max_layer_height",0.);
    }
    const double filament_diameter = first(config,"filament_diameter",1.75), flow_ratio = first(config,"filament_flow_ratio",1.);
    if (!std::isfinite(filament_diameter) || !std::isfinite(flow_ratio) || filament_diameter <= 0. || flow_ratio <= 0.)
        invalid(_u8L("Filament diameter and flow ratio must be positive and finite"));
    s.min_deposited_volume = .25*std::acos(-1.)*filament_diameter*filament_diameter*1e-6/flow_ratio;
    s.width = config.get_abs_value("line_width",s.nozzle_diameter);
    if (s.width == 0.) s.width = s.nozzle_diameter*1.125;
    // All roles use one explicit projected width. Role-specific widths would
    // change the erosion, physical cross-section, and spacing together.
    for (auto key : {"outer_wall_line_width","inner_wall_line_width","sparse_infill_line_width", "internal_solid_infill_line_width","top_surface_line_width", "initial_layer_line_width"}) {
        const double width = config.get_abs_value(key,s.nozzle_diameter);
        if (width != 0. && std::abs(width-s.width) > 1e-5)
            invalid(Slic3r::format(_u8L("%1% must equal line_width (or be zero for automatic width)"),key));
    }
    s.walls = config.opt_int("wall_loops");
    s.top_layers = config.opt_int("top_shell_layers");
    s.bottom_layers = config.opt_int("bottom_shell_layers");
    s.density = scalar(config,"sparse_infill_density",15.)*.01;
    s.infill_angle = scalar(config,"infill_direction",45.);
    s.max_z = scalar(config,"printable_height",250.);
    s.clearance_radius = scalar(config,"extruder_clearance_radius",0.);
    s.clearance_height = std::min({scalar(config,"extruder_clearance_height_to_rod",0.),
                                  scalar(config,"extruder_clearance_height_to_lid",0.), scalar(config,"nozzle_height",0.)});
    if (s.clearance_radius <= 0 || s.clearance_height <= 0) invalid(_u8L("positive toolhead clearance radius and height are required"));
    std::vector<Object> objects;
    const Vec3d origin = print.get_plate_origin();
    Polygon bed_polygon = Polygon::new_scale(config.option<ConfigOptionPoints>("printable_area")->values);
    bed_polygon.make_counter_clockwise();
    if (!polygon_is_convex(bed_polygon))
        invalid(_u8L("nonplanar travel currently requires a convex printable area"));
    for (const auto *po : print.objects()) {
        const auto *mo = po->model_object();
        for (const auto &key : mo->config.get().keys())
            if (key != "extruder" || mo->config.get().opt_int(key) > 1) invalid(_u8L("per-object overrides are unsupported"));
        if (mo->volumes.size() != 1) invalid(_u8L("use one watertight model volume per object"));
        for (const auto &instance : po->instances()) {
            Object object;
            object.id = po->id().id;
            for (const auto *volume : mo->volumes) {
                if (!volume->is_model_part()) invalid(_u8L("modifier and negative volumes are unsupported; use a watertight mesh with native cavities"));
                for (const auto &key : volume->config.get().keys())
                    if (key != "extruder" || volume->config.get().opt_int(key) > 1) invalid(_u8L("per-volume overrides are unsupported"));
                const Transform3d transform = instance.model_instance->get_matrix()*volume->get_matrix();
                const auto &mesh = volume->mesh().its;
                validate_mesh(mesh);
                for (const auto &face : mesh.indices) {
                    Triangle triangle;
                    for (int n = 0; n < 3; ++n) {
                        const Vec3d p = transform*mesh.vertices[face[n]].cast<double>()-origin;
                        if (!bed_polygon.contains(Point(scale_(p.x()),scale_(p.y()))))
                            invalid(_u8L("model is outside the printable area"));
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
    parser.set("current_extruder",0);
    parser.set("total_layer_count",int(job.layers.size()));
    size_t last_printed = 0;
    const Path *finishing_path = nullptr;
    double max_printed_z = 0.;
    for (size_t k = 0; k < job.layers.size(); ++k)
        for (const Path &path : job.layers[k].paths)
            if (std::any_of(path.segments.begin(),path.segments.end(),[](const Segment &s) { return s.volume > 0.; })) {
                last_printed = k; finishing_path = &path;
                max_printed_z = std::max(max_printed_z,job.layers[k].max_z);
            }
    parser.set("layer_num",0);
    parser.set("layer_z",job.layers.front().max_z);
    parser.set("max_layer_z",max_printed_z);
    parser.set("max_print_z",std::ceil(max_printed_z));
    parser.set("max_print_height",job.settings.max_z);
    parser.set("z_offset",0.);
    parser.set("model_name",print.get_model_name());
    parser.set("plate_number",print.get_plate_number_formatted());
    BoundingBoxf first_layer_box;
    for (const auto &path : job.layers.front().paths) {
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
        << "; nonplanar_bead_model = vertical_band_projected_width_v1\n; FEATURE: Custom\n";
    out << writer.set_bed_temperature(initial_bed->values.front(),false)
        << writer.set_temperature(initial_nozzle->values.front(),false)
        << startup << '\n'
        << writer.set_bed_temperature(initial_bed->values.front(),true)
        << writer.set_temperature(initial_nozzle->values.front(),true)
        << filament_start << '\n' << writer.preamble() << writer.reset_e(true);
    double printed_max = job.stack.bounds(job.stack.bed).second;
    double volume = 0.;
    auto clear_travel = [&](const Vec3d &target) {
        const double safe_z = std::max(printed_max,target.z()) + std::max(1.,job.stack.limits.max_gap)+tol;
        if (safe_z > job.settings.max_z) invalid(_u8L("collision-free travel exceeds machine Z height"));
        if (target.z()+job.settings.clearance_height <= printed_max+tol)
            invalid(_u8L("toolhead clearance cannot reach this path below already printed geometry"));
        out << writer.retract() << writer.set_speed(xy_travel*60.)
            << writer.travel_to_z(safe_z,"nonplanar clearance",true)
            << writer.travel_to_xyz({target.x(),target.y(),safe_z},"nonplanar travel",true)
            << writer.travel_to_xyz(target,"nonplanar approach",true) << writer.unretract();
    };
    for (size_t k = 0; k < job.layers.size(); ++k) {
        if (print.canceled()) throw CanceledException();
        const Layer &layer = job.layers[k];
        if (k == 1) out << writer.set_bed_temperature(bed->values.front(),false) << writer.set_temperature(nozzle->values.front(),false);
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
                if (end.z()+job.settings.clearance_height <= printed_max+tol)
                    invalid(_u8L("toolhead clearance fails along a nonplanar deposition path"));
                // Native metadata drives modeled preview dimensions per segment.
                out << "; LINE_WIDTH: " << segment.physical_width << "\n; LAYER_HEIGHT: " << segment.normal_height
                    << "\n;NONPLANAR_BAND:" << job.settings.width << ',' << segment.vertical_height << ','
                    << segment.normal_height << ',' << segment.cross_slope << '\n';
                out << writer.set_speed(feed*60.) << writer.extrude_to_xyz(end,segment.volume*ratio/filament_area);
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
