#pragma once

#include "LayerStack.hpp"
#include <array>
#include <cstdint>
#include <functional>
#include <vector>

namespace Slic3r::Nonplanar {

struct Triangle { std::array<XYZ,3> vertices; };
struct Object { std::vector<Triangle> triangles; uint64_t id = 0; };
enum class Role { OuterWall, InnerWall, SparseInfill, SolidInfill, TopSkin, BottomSkin };

// Dimensions describe the actual swept band, not a planar bead relocated in Z.
struct Segment {
    XYZ end;
    double volume, vertical_height, normal_height, physical_width, cross_slope, max_top_z;
};
struct Path {
    XYZ start;
    Role role;
    std::vector<Segment> segments;
    uint64_t object_id = 0;
    double projected_area = 0.;
};
struct Layer {
    std::vector<Path> paths;
    double min_z = 0., max_z = 0.;
};
struct Settings {
    int walls = 2, top_layers = 3, bottom_layers = 3;
    double width = 0.45, density = 0.15;
    double infill_angle = 0.;
    double nozzle_diameter = 0.4, max_z = 250.;
    double min_normal_height = 0., max_normal_height = 0.;
    double min_deposited_volume = 0.;
    double clearance_radius = 0., clearance_height = 0.;
};

struct Job {
    LayerStack stack;
    Settings settings;
    std::vector<Layer> layers;
};

// Native original-mesh/spline contours, 2D toolpaths and local gap integration.
// Cancellation is supplied by Print, without a GUI or service dependency.
std::vector<Layer> slice(const LayerStack &, const std::vector<Object> &, const Settings &,
                         const std::function<void()> &cancel);

} // namespace Slic3r::Nonplanar
