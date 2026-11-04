#pragma once

#include "Slicer.hpp"
#include <memory>
#include <string>

namespace Slic3r {
class Model;
class Print;
struct GCodeProcessorResult;
namespace Nonplanar {
// Model metadata is the durable mode switch. A plain model has no nonplanar job.
std::string source(const Model &, int plate_index);
void attach(Model &, const LayerStack &);
void detach(Model &);
std::shared_ptr<Job> prepare(Print &, const std::string &source);
std::string export_gcode(Print &, const Job &, const std::string &path, GCodeProcessorResult *);
} // namespace Nonplanar
} // namespace Slic3r
