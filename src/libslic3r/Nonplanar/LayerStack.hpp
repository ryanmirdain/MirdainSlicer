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

// In one shared grid cell: c + ax*u + ay*v + axy*u*v, in millimeters.
// u and v are physical offsets from the cell's lower left corner.
struct Patch {
    double x, y, c, ax, ay, axy;
    double value(XY p) const;
    XY gradient(XY p) const;
};

struct Surface {
    bool plane = true;
    double a = 0., b = 0., c = 0.;
    std::vector<double> z; // row-major, shared grid
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
