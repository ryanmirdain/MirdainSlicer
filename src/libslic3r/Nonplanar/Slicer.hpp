#pragma once

#include "LayerStack.hpp"
#include "../Point.hpp"
#include "../Polygon.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>
#include <utility>

namespace Slic3r::Nonplanar {
struct Triangle {
    std::array<XYZ, 3> vertices;
    std::array<size_t, 3> vertex_ids;
};
struct Object { std::vector<Triangle> triangles; };
struct Job {
    explicit Job(std::shared_ptr<const LayerStack> layers) : m_stack(std::move(layers)) {}
    const LayerStack &stack() const { return *m_stack; }
    double first_layer_height, layer_height;
    double layer_z(size_t index) const { return first_layer_height + index * layer_height; }
    Surface surface_at(double slicer_z) const;
    // Present when the sampling surface is a vertical copy of its layer top.
    std::optional<double> surface_offset(double slicer_z) const;
    Vec3d to_model(const Vec3d &slicer_point) const;
private:
    std::shared_ptr<const LayerStack> m_stack;
};

// One immutable mesh index per placed volume. Only finished contours are
// cached; tessellation storage lives inside each parallel layer task.
class MeshSlicer {
public:
    MeshSlicer(const Job &, const Object &, const std::function<void()> &cancel);
    ~MeshSlicer();
    std::vector<Polygons> cross_sections(const std::vector<float> &zs, const std::function<void()> &cancel,
                                       double reference_offset, bool reference);
    double index_seconds() const;
    double contour_seconds() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
} // namespace Slic3r::Nonplanar
