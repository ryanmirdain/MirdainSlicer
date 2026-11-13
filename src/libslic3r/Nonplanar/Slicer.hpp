#pragma once

#include "LayerStack.hpp"
#include "../Point.hpp"
#include "../Polygon.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Slic3r::Nonplanar {
struct Triangle {
    std::array<XYZ, 3> vertices;
    std::array<size_t, 3> vertex_ids;
};
struct Object { std::vector<Triangle> triangles; };
struct Job {
    LayerStack stack;
    double first_layer_height, layer_height;
    double layer_z(size_t index) const { return first_layer_height + index * layer_height; }
    Surface surface_at(double slicer_z) const;
    Vec3d to_model(const Vec3d &slicer_point) const;
};

// Intersect unchanged triangles at the inverse images of Orca's slicing Zs.
std::vector<Polygons> cross_sections(const Job &, const Object &, const std::vector<float> &zs,
                                    const std::function<void()> &cancel);
} // namespace Slic3r::Nonplanar
