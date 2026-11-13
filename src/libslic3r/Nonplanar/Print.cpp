#include "Print.hpp"
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
    job->first_layer_height = print.config().initial_layer_print_height.value;
    job->layer_height = print.default_object_config().layer_height.value;
    return job;
}

std::vector<ExPolygons> slice_volume(const indexed_triangle_set &mesh, const std::vector<float> &zs,
                                    const MeshSlicingParamsEx &params, const SliceContext &context,
                                    const std::function<void()> &cancel)
{
    Object original;
    original.triangles.reserve(mesh.indices.size());
    for (const auto &face : mesh.indices) {
        Triangle triangle;
        for (size_t corner = 0; corner < 3; ++corner) {
            Vec3d p = params.trafo * mesh.vertices[face[corner]].cast<double>();
            p.head<2>() += context.xy_offset;
            triangle.vertices[corner] = {p.x(), p.y(), p.z()};
            triangle.vertex_ids[corner] = size_t(face[corner]);
        }
        original.triangles.push_back(triangle);
    }
    auto sections = cross_sections(context.job, original, zs, cancel);
    for (size_t layer = 0; layer < sections.size(); ++layer) {
        const auto mode = layer < params.slicing_mode_normal_below_layer ? params.mode_below : params.mode;
        for (auto &polygon : sections[layer]) {
            polygon.translate(Point(-scaled<coord_t>(context.xy_offset.x()), -scaled<coord_t>(context.xy_offset.y())));
            if (mode == MeshSlicingParams::SlicingMode::Positive || mode == MeshSlicingParams::SlicingMode::PositiveLargestContour)
                polygon.make_counter_clockwise();
        }
    }
    return slice_mesh_ex(std::move(sections), params, cancel);
}
} // namespace Slic3r::Nonplanar
