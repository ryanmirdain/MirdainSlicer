#pragma once

#include "LayerStack.hpp"
#include <array>
#include <functional>
#include <vector>

namespace Slic3r::Nonplanar {

struct Triangle { std::array<XYZ,3> vertices; };
struct Object { std::vector<Triangle> triangles; };
enum class Role { OuterWall, InnerWall, SparseInfill, SolidInfill, TopSkin, BottomSkin };

// Dimensions describe the actual swept band, not a planar bead relocated in Z.
struct Segment {
    XYZ end;
    double volume, vertical_height, normal_height, physical_width, cross_slope;
};
struct Path {
    XYZ start;
    Role role;
    std::vector<Segment> segments;
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
    double clearance_radius = 0., clearance_height = 0.;
};

struct Job {
    LayerStack stack;
    Settings settings;
    std::vector<Layer> layers;
};

// Native original-mesh/bilinear intersection, band booleans and swept-volume paths.
// Cancellation is supplied by Print, without a GUI or service dependency.
std::vector<Layer> slice(const LayerStack &, const std::vector<Object> &, const Settings &,
                         const std::function<void()> &cancel);

} // namespace Slic3r::Nonplanar
