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
class ModelVolume;
class Print;
struct GCodeProcessorResult;
namespace Nonplanar {
std::shared_ptr<const LayerStack> layer_stack(const Model &, int plate_index,
                                              const std::shared_ptr<const LayerStack> &previous = {});
void attach(Model &, const std::shared_ptr<const LayerStack> &);
void detach(Model &);
std::shared_ptr<Job> prepare(Print &, const std::shared_ptr<const LayerStack> &);
struct SliceCache;
std::shared_ptr<SliceCache> make_slice_cache();
struct SliceContext {
    const Job &job;
    Vec2d xy_offset;
    bool reference = false;
    double reference_offset = 0.;
    std::shared_ptr<SliceCache> cache = make_slice_cache();
};
std::vector<ExPolygons> slice_volume(const ModelVolume &, const std::vector<float> &zs,
                                    const MeshSlicingParamsEx &, const SliceContext &,
                                    const std::function<void()> &cancel);
void log_slice_timings(const SliceContext &);
// Historical inverse: map model movements, retaining total E and feed rates.
void unwarp_gcode(Print &, const Job &, const std::string &path, GCodeProcessorResult *);
} // namespace Nonplanar
} // namespace Slic3r
