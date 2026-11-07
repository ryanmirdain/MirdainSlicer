#include "Print.hpp"
#include "../GCodeReader.hpp"
#include "../GCode/GCodeProcessor.hpp"
#include "../Exception.hpp"
#include "../I18N.hpp"
#include "../LocalesUtils.hpp"
#include "../Print.hpp"
#include "../Point.hpp"
#include "../Utils.hpp"
#include <boost/nowide/fstream.hpp>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>
#include <utility>

namespace Slic3r::Nonplanar {
namespace {
constexpr const char *inverse_marker = "; Nonplanar inverse: absolute XYZ, relative E";

std::string trimmed_upper(std::string value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    value = value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return char(std::toupper(c)); });
    return value;
}

using Parameters = std::map<char, double>;
Parameters parameters(const GCodeReader::GCodeLine &line)
{
    Parameters values;
    for (char axis = 'A'; axis <= 'Z'; ++axis) {
        const auto word = line.axis_pos(axis);
        if (word.empty()) continue;
        char *end;
        const double value = std::strtod(word.data() + 1, &end);
        if (end != word.data() + 1) values[axis] = value;
    }
    return values;
}

std::string extras(const GCodeReader::GCodeLine &line)
{
    const auto command = line.cmd();
    const auto offset = size_t(command.data() - line.raw().data()) + command.size();
    const auto comment = line.raw().find(';', offset);
    std::istringstream words(line.raw().substr(offset, comment - offset));
    std::string word, result;
    while (words >> word)
        if (std::string("XYZEIJR").find(char(std::toupper(static_cast<unsigned char>(word.front())))) == std::string::npos)
            result += " " + word;
    return result;
}

std::vector<Vec3d> motion_points(const Vec3d &from, const Vec3d &to, const Parameters &p,
                                const std::string &command, bool absolute_center, double maximum_length)
{
    if (command == "G0" || command == "G1") {
        const size_t count = std::max(size_t(1), size_t(std::ceil((to - from).norm() / maximum_length)));
        std::vector<Vec3d> points;
        points.reserve(count);
        for (size_t n = 1; n <= count; ++n)
            points.push_back(from + double(n) / double(count) * (to - from));
        return points;
    }
    const bool clockwise = command == "G2";
    std::vector<Vec2d> centers;
    if (p.count('R')) {
        const Vec2d delta = to.head<2>() - from.head<2>();
        const double length = delta.norm(), radius = std::abs(p.at('R'));
        if (length == 0. || length > 2. * radius + 1e-6)
            throw std::invalid_argument(_u8L("Invalid radius-defined G-code arc"));
        const Vec2d midpoint = .5 * (from.head<2>() + to.head<2>());
        const Vec2d perpendicular(-delta.y() / length, delta.x() / length);
        const double distance = std::sqrt(std::max(0., radius * radius - length * length * .25));
        centers = {midpoint + distance * perpendicular, midpoint - distance * perpendicular};
    } else if (p.count('I') || p.count('J')) {
        Vec2d center(p.count('I') ? p.at('I') : absolute_center ? from.x() : 0.,
                     p.count('J') ? p.at('J') : absolute_center ? from.y() : 0.);
        if (!absolute_center) center += from.head<2>();
        centers.push_back(center);
    } else {
        throw std::invalid_argument(_u8L("An XY arc requires I/J or R"));
    }
    double angle = 0., sweep = 0., radius = 0.;
    Vec2d center;
    const double pi = std::acos(-1.);
    for (const auto &candidate : centers) {
        center = candidate;
        radius = (from.head<2>() - center).norm();
        angle = std::atan2(from.y() - center.y(), from.x() - center.x());
        sweep = std::atan2(to.y() - center.y(), to.x() - center.x()) - angle;
        if (clockwise && sweep >= 0.) sweep -= 2. * pi;
        if (!clockwise && sweep <= 0.) sweep += 2. * pi;
        if (!p.count('R') || (std::abs(sweep) > pi + 1e-8) == (p.at('R') < 0.)) break;
    }
    const double length = std::hypot(radius * std::abs(sweep), to.z() - from.z());
    const size_t count = std::max(size_t(1), size_t(std::ceil(length / maximum_length)));
    std::vector<Vec3d> points;
    points.reserve(count);
    for (size_t n = 1; n <= count; ++n) {
        const double fraction = double(n) / double(count), theta = angle + sweep * fraction;
        points.emplace_back(center.x() + radius * std::cos(theta), center.y() + radius * std::sin(theta),
                            from.z() + fraction * (to.z() - from.z()));
    }
    points.back() = to;
    return points;
}
} // namespace

void unwarp_gcode(Print &print, const Job &job, const std::string &path, GCodeProcessorResult *result)
{
    CNumericLocalesSetter locale;
    std::vector<std::string> lines;
    boost::nowide::ifstream input(path, std::ios::binary);
    input.exceptions(std::ios::badbit);
    std::string text;
    while (std::getline(input, text)) lines.push_back(text + "\n");
    input.close();

    size_t first_part = lines.size(), first = lines.size(), last = lines.size();
    bool already_mapped = false;
    for (size_t index = 0; index < lines.size(); ++index) {
        const auto line = trimmed_upper(lines[index]);
        if (line.empty() || line.front() != ';') continue;
        const auto tag = trimmed_upper(line.substr(1));
        if (tag == trimmed_upper(std::string(inverse_marker).substr(1))) already_mapped = true;
        if (first_part == lines.size() && (tag.rfind("START PRINTING OBJECT,", 0) == 0 || tag.rfind("PRINTING OBJECT ", 0) == 0))
            first_part = index;
        if (index > first_part && first == lines.size() && (tag.rfind("FEATURE:", 0) == 0 || tag.rfind("TYPE:", 0) == 0))
            first = index;
        if (index > first && tag == "MACHINE_END_GCODE_START") { last = index; break; }
    }
    if (!already_mapped) {
        if (first == lines.size() || last == lines.size())
            throw std::invalid_argument(_u8L("Orca part/feature and machine-end markers are required for nonplanar inverse mapping"));
        const std::string temporary = path + ".unwarp";
        boost::nowide::ofstream output(temporary, std::ios::binary);
        output.exceptions(std::ios::failbit | std::ios::badbit);
        output.imbue(std::locale::classic());
        output << std::fixed << std::setprecision(6);
        Vec3d position = Vec3d::Constant(std::numeric_limits<double>::quiet_NaN());
        double e_position = 0.;
        bool absolute_xyz = true, absolute_e = true, absolute_center = false;
        std::string plane = "G17";
        Vec2d extruder_offset = print.config().extruder_offset.get_at(0);
        GCodeReader reader;
        auto emit_point = [&](const Vec3d &point, bool z_only = false) {
            // Orca subtracts the extruder offset from commanded XY. Surface
            // lookup needs the nozzle's model position, not carriage position.
            Vec3d model_point = point;
            model_point.head<2>() += extruder_offset;
            const Vec3d mapped = job.to_model(model_point);
            if (!z_only) output << " X" << point.x() << " Y" << point.y();
            output << " Z" << std::max(.01, mapped.z());
        };
        for (size_t index = 0; index < last; ++index) {
            if ((index & 1023) == 0 && print.canceled()) throw CanceledException();
            const bool body = index >= first;
            if (index == first) {
                if (!position.allFinite()) throw std::invalid_argument(_u8L("Model motion starts before all XYZ coordinates are established"));
                output << inverse_marker << "\nG90\nM83\nG1";
                emit_point(position, true);
                output << '\n';
            }
            reader.parse_line(lines[index], [&](GCodeReader &, const GCodeReader::GCodeLine &line) {
                const std::string command = trimmed_upper(std::string(line.cmd()));
                const auto p = parameters(line);
                if (command.size() > 1 && command.front() == 'T' && std::isdigit(static_cast<unsigned char>(command[1]))) {
                    const size_t filament = std::stoul(command.substr(1));
                    if (filament < print.config().filament_map.size())
                        extruder_offset = print.config().extruder_offset.get_at(print.config().filament_map.get_at(filament) - 1);
                }
                if (command == "G90" || command == "G91" || command == "M82" || command == "M83") {
                    if (command[0] == 'G') absolute_xyz = command == "G90";
                    else absolute_e = command == "M82";
                    if (!body) output << lines[index];
                    return;
                }
                if (command == "G90.1" || command == "G91.1") absolute_center = command == "G90.1";
                if (command == "G17" || command == "G18" || command == "G19") plane = command;
                if (command == "G20") throw std::invalid_argument(_u8L("Nonplanar inverse mapping requires millimeter G-code"));
                if (body && (command == "G53" || command == "G54" || command == "G55" || command == "G56" ||
                             command == "G57" || command == "G58" || command == "G59" || command == "G28"))
                    throw std::invalid_argument(_u8L("Work-coordinate changes and homing inside model movements are unsupported"));
                if (command == "G28") {
                    bool any = false;
                    for (char axis : std::string("XYZ")) any |= line.has(axis);
                    for (size_t axis = 0; axis < 3; ++axis)
                        if (!any || line.has("XYZ"[axis])) position[axis] = std::numeric_limits<double>::quiet_NaN();
                }
                if (command == "G92") {
                    for (size_t axis = 0; axis < 3; ++axis) if (p.count("XYZ"[axis])) {
                        if (body) throw std::invalid_argument(_u8L("G92 XYZ inside model movements is unsupported"));
                        position[axis] = p.at("XYZ"[axis]);
                    }
                    if (p.count('E')) e_position = p.at('E');
                }
                if (command != "G0" && command != "G1" && command != "G2" && command != "G3") {
                    output << lines[index]; return;
                }
                const Vec3d from = position;
                for (size_t axis = 0; axis < 3; ++axis)
                    if (p.count("XYZ"[axis])) position[axis] = (absolute_xyz ? 0. : position[axis]) + p.at("XYZ"[axis]);
                const bool has_e = p.count('E');
                const double delta_e = has_e ? p.at('E') - (absolute_e ? e_position : 0.) : 0.;
                e_position += delta_e;
                if (!body) { output << lines[index]; return; }
                const bool has_xyz = p.count('X') || p.count('Y') || p.count('Z');
                const std::string other = extras(line);
                const std::string comment = line.comment().empty() ? "" : ";" + std::string(line.comment());
                if (!has_xyz && (command == "G0" || command == "G1")) {
                    output << command;
                    if (has_e) output << " E" << delta_e;
                    output << other;
                    if (!comment.empty()) output << " " << comment;
                    output << '\n'; return;
                }
                if (!from.allFinite() || !position.allFinite())
                    throw std::invalid_argument(_u8L("Model motion requires known XYZ coordinates"));
                if ((command == "G2" || command == "G3") && plane != "G17")
                    throw std::invalid_argument(_u8L("Only XY-plane arcs can be inverse mapped"));
                const auto points = motion_points(from, position, p, command, absolute_center, job.stack.limits.max_segment_length);
                double allocated = 0.;
                for (size_t n = 0; n < points.size(); ++n) {
                    output << (command == "G0" ? "G0" : "G1");
                    emit_point(points[n]);
                    if (has_e) {
                        const double target = std::round(delta_e * double(n + 1) / double(points.size()) * 1e6) / 1e6;
                        output << " E" << target - allocated;
                        allocated = target;
                    }
                    if (n == 0) output << other;
                    if (n + 1 == points.size() && !comment.empty()) output << " " << comment;
                    output << '\n';
                }
            });
        }
        output << "; Restore Orca modes for shutdown\nG92 E" << e_position << '\n'
               << (absolute_xyz ? "G90\n" : "G91\n") << (absolute_e ? "M82\n" : "M83\n");
        for (size_t index = last; index < lines.size(); ++index) output << lines[index];
        output.close();
        if (const auto error = rename_file(temporary, path)) throw std::system_error(error, "Replacing inverse-mapped G-code");
    }
    GCodeProcessor processor;
    processor.apply_config(print.config());
    processor.set_print(&print);
    const Vec3d origin = print.get_plate_origin();
    processor.set_xy_offset(origin.x(), origin.y());
    processor.process_file(path, [&print] { if (print.canceled()) throw CanceledException(); });
    print.print_statistics().estimated_normal_print_time = processor.get_time_dhm(PrintEstimatedStatistics::ETimeMode::Normal);
    if (result) {
        auto &&mapped = processor.extract_result();
        result->moves = std::move(mapped.moves);
        result->lines_ends = std::move(mapped.lines_ends);
        result->print_statistics.modes = std::move(mapped.print_statistics.modes);
        result->print_statistics.total_travel_distance = mapped.print_statistics.total_travel_distance;
        result->print_statistics.total_travel_moves = mapped.print_statistics.total_travel_moves;
        result->skippable_part_time = std::move(mapped.skippable_part_time);
    }
}
} // namespace Slic3r::Nonplanar
