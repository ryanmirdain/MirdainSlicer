// Copyright (c) 2026. Distributed under the GNU AGPLv3, like OrcaSlicer.
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Nonplanar {

struct XY { double x, y; };
struct XYZ { double x, y, z; };

// One cubic Bezier span, or two planar triangles when linear is set.
struct Patch {
    double x, y, dx, dy;
    std::array<double,16> control {};
    bool linear = false; // Two planar triangles sharing the cell's lower-left/upper-right diagonal.
    double value(XY p) const;
    XY gradient(XY p) const;
    std::pair<double,double> bounds() const;
    double linear_error() const;
    std::array<double,7> line(XY a, XY b) const;
};

struct Surface {
    bool plane = true;
    double a = 0., b = 0., c = 0.;
    std::vector<double> z; // row-major input samples
    std::vector<double> gx, gy; // shared C1 material-constrained derivatives
    double level = 0.;
    std::vector<Patch> patches; // positive-weight, shared-knot spline spans
    double min_z = 0., max_z = 0.; // derived grid bounds; input surfaces are immutable
};

struct Limits {
    double geometry_tolerance, chord_tolerance, max_segment_length;
    double min_gap, max_gap, max_slope_degrees;
};

class LayerStack {
public:
    static LayerStack parse(const std::string &json);
    static LayerStack read(const std::string &filename);

    int plate_index = 0;
    std::vector<double> x, y;
    Surface bed;
    std::vector<Surface> layers;
    std::vector<XYZ> field_vertices;
    std::vector<double> field_values;
    std::vector<std::array<size_t,3>> field_triangles;
    Limits limits {};
    std::string source; // validated, canonical JSON; embedded verbatim in projects

    bool covers(XY p) const;
    std::pair<size_t, size_t> cell(XY p) const;
    Patch patch(const Surface &, size_t i, size_t j) const;
    double value(const Surface &, XY p) const;
    XY gradient(const Surface &, XY p) const;
    std::pair<double, double> bounds(const Surface &) const;
    double vertical_gap(size_t layer, XY p) const;
    double normal_spacing(size_t layer, XY p) const;
};

inline constexpr const char *metadata_key = "orca_nonplanar_layers_v1";
inline constexpr const char *version_key = "orca_nonplanar_schema_version";

} // namespace Slic3r::Nonplanar
