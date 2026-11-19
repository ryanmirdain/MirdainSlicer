#include "Print.hpp"
#include "../ClipperUtils.hpp"
#include "../I18N.hpp"
#include "../format.hpp"
#include "../Print.hpp"
#include "../Model.hpp"
#include "../PrintConfig.hpp"
#include "../TriangleMesh.hpp"
#include "../Polygon.hpp"
#include "../Point.hpp"
#include "../Exception.hpp"
#include <Eigen/Geometry>
#include <boost/log/trivial.hpp>
#include <chrono>
#include <map>
#include <mutex>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>
#include <cmath>
#include <memory>
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

} // namespace

std::shared_ptr<const LayerStack> layer_stack(const Model &model, int plate,
                                              const std::shared_ptr<const LayerStack> &previous)
{
    if (!model.model_info) return {};
    const auto &items = model.model_info->metadata_items;
    const auto it = items.find(metadata_key);
    if (it == items.end()) return {};
    const auto version = items.find(version_key);
    if (version == items.end() || (version->second != "1" && version->second != "2")) invalid(_u8L("unsupported or missing embedded nonplanar schema version"));
    auto stack = previous && previous->source == it->second ? previous : model.model_info->nonplanar_layers;
    if (!stack || stack->source != it->second)
        stack = std::make_shared<LayerStack>(LayerStack::parse(it->second));
    return stack->plate_index == plate ? stack : nullptr;
}

void attach(Model &model, const std::shared_ptr<const LayerStack> &stack)
{
    model.model_info = model.model_info ? std::make_shared<ModelInfo>(*model.model_info) : std::make_shared<ModelInfo>();
    auto &items = model.model_info->metadata_items;
    items[metadata_key] = stack->source;
    items[version_key] = std::to_string(stack->version);
    model.model_info->nonplanar_layers = stack;
    items[placement_key] = placement(model);
}
void detach(Model &model)
{
    if (!model.model_info) return;
    model.model_info = std::make_shared<ModelInfo>(*model.model_info);
    for (auto key : {metadata_key,version_key,placement_key}) model.model_info->metadata_items.erase(key);
    model.model_info->nonplanar_layers.reset();
}

std::shared_ptr<Job> prepare(Print &print, const std::shared_ptr<const LayerStack> &stack)
{
    auto job = std::make_shared<Job>(stack);
    job->first_layer_height = print.config().initial_layer_print_height.value;
    job->layer_height = print.default_object_config().layer_height.value;
    return job;
}

struct SliceCache {
    std::mutex mutex;
    std::map<const ModelVolume *,std::shared_ptr<MeshSlicer>> meshes;
    double transform_seconds = 0., assembly_seconds = 0., mask_seconds = 0.;
};
std::shared_ptr<SliceCache> make_slice_cache() { return std::make_shared<SliceCache>(); }

std::vector<ExPolygons> slice_volume(const ModelVolume &volume, const std::vector<float> &zs,
                                    const MeshSlicingParamsEx &params, const SliceContext &context,
                                    const std::function<void()> &cancel)
{
    using Clock = std::chrono::steady_clock;
    std::shared_ptr<MeshSlicer> slicer;
    {
        std::lock_guard<std::mutex> lock(context.cache->mutex);
        auto &cached = context.cache->meshes[&volume];
        if (!cached) {
            const auto started = Clock::now();
            Object original;
            const auto &mesh = volume.mesh().its;
            const bool flipped = params.trafo.rotation().determinant() < 0.;
            original.triangles.reserve(mesh.indices.size());
            for (auto face : mesh.indices) {
                if (flipped) std::swap(face[1],face[2]);
                Triangle triangle;
                for (size_t corner = 0; corner < 3; ++corner) {
                    Vec3d p = params.trafo * mesh.vertices[face[corner]].cast<double>();
                    p.head<2>() += context.xy_offset;
                    triangle.vertices[corner] = {p.x(),p.y(),p.z()};
                    triangle.vertex_ids[corner] = size_t(face[corner]);
                }
                original.triangles.push_back(triangle);
            }
            context.cache->transform_seconds += std::chrono::duration<double>(Clock::now()-started).count();
            cached = std::make_shared<MeshSlicer>(context.job,original,cancel);
        }
        slicer = cached;
    }
    auto sections = slicer->cross_sections(zs,cancel,context.reference_offset,context.reference);
    const auto assembly_started = Clock::now();
    for (size_t layer = 0; layer < sections.size(); ++layer) {
        const auto mode = layer < params.slicing_mode_normal_below_layer ? params.mode_below : params.mode;
        for (auto &polygon : sections[layer]) {
            polygon.translate(Point(-scaled<coord_t>(context.xy_offset.x()), -scaled<coord_t>(context.xy_offset.y())));
            if (mode == MeshSlicingParams::SlicingMode::Positive || mode == MeshSlicingParams::SlicingMode::PositiveLargestContour)
                polygon.make_counter_clockwise();
        }
    }
    auto result = slice_mesh_ex(std::move(sections), params, cancel);
    const auto mask_started = Clock::now();
    if (!context.reference && context.job.stack().masked) {
        const Point shift(scaled<coord_t>(context.xy_offset.x()),scaled<coord_t>(context.xy_offset.y()));
        for (size_t n = 0; n < result.size(); ++n) {
            const size_t k = std::min(context.job.stack().layers.size()-1,
                size_t(std::max(0.,std::ceil((double(zs[n])-context.job.first_layer_height)/context.job.layer_height))));
            if (!context.job.stack().has_clip_mask[k]) continue;
            ExPolygons mask = context.job.stack().clip_masks[k];
            for (auto &polygon : mask) polygon.translate(-shift);
            result[n] = intersection_ex(result[n],mask);
        }
    }
    const auto finished = Clock::now();
    {
        std::lock_guard<std::mutex> lock(context.cache->mutex);
        context.cache->assembly_seconds += std::chrono::duration<double>(mask_started-assembly_started).count();
        context.cache->mask_seconds += std::chrono::duration<double>(finished-mask_started).count();
    }
    return result;
}

void log_slice_timings(const SliceContext &context)
{
    std::lock_guard<std::mutex> lock(context.cache->mutex);
    double index = 0., contours = 0.;
    for (const auto &[volume,mesh] : context.cache->meshes) {
        index += mesh->index_seconds();
        contours += mesh->contour_seconds();
    }
    // Per-volume wall times; independent volumes may overlap in the TBB pool.
    BOOST_LOG_TRIVIAL(info) << "Nonplanar timings: " << nlohmann::json({
        {"mesh_count",context.cache->meshes.size()},
        {"transform_seconds",context.cache->transform_seconds},
        {"index_seconds",index},{"contour_seconds",contours},
        {"assembly_seconds",context.cache->assembly_seconds},
        {"mask_seconds",context.cache->mask_seconds}}).dump();
}
} // namespace Slic3r::Nonplanar
