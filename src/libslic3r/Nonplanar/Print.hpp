#pragma once

#include "Slicer.hpp"
#include "../ExPolygon.hpp"
#include "../TriangleMesh.hpp"
#include "../TriangleMeshSlicer.hpp"
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Slic3r {
class Model;
class Print;
struct GCodeProcessorResult;
namespace Nonplanar {
std::string source(const Model &, int plate_index);
void attach(Model &, const LayerStack &);
void detach(Model &);
std::shared_ptr<Job> prepare(Print &, const std::string &source);
struct SliceContext {
    const Job &job;
    Vec2d xy_offset;
};
std::vector<ExPolygons> slice_volume(const indexed_triangle_set &, const std::vector<float> &zs,
                                    const MeshSlicingParamsEx &, const SliceContext &,
                                    const std::function<void()> &cancel);
// Historical inverse: map model movements, retaining total E and feed rates.
void unwarp_gcode(Print &, const Job &, const std::string &path, GCodeProcessorResult *);
} // namespace Nonplanar
} // namespace Slic3r
