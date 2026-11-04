#include "LayerStack.hpp"

#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace Slic3r::Nonplanar {
namespace {
using Json = nlohmann::json;
[[noreturn]] void invalid(const std::string &message)
{ throw std::invalid_argument("Nonplanar layers: " + message); }

double number(const Json &j, const std::string &label)
{
    if (!j.is_number() || !std::isfinite(j.get<double>())) invalid(label + " must be finite");
    return j.get<double>();
}

std::vector<double> axis(const Json &j, const std::string &name)
{
    if (!j.is_array() || j.size() < 2) invalid("grid." + name + " needs at least two values");
    std::vector<double> out;
    for (const auto &v : j) {
        const double n = number(v, "grid." + name);
        // Clipper coordinates use 1e-6 mm. Larger coordinates are not printer coordinates.
        if (std::abs(n) > 1e7) invalid("grid coordinates exceed the supported physical range");
        if (!out.empty() && n <= out.back()) invalid("grid." + name + " must be strictly increasing");
        out.push_back(n);
    }
    return out;
}

Surface surface(const Json &j, size_t nx, size_t ny, const std::string &label)
{
    Surface s;
    const std::string type = j.at("type").get<std::string>();
    if (type == "plane") {
        s.a = number(j.at("a"), label + ".a");
        s.b = number(j.at("b"), label + ".b");
        s.c = number(j.at("c"), label + ".c");
    } else if (type == "height_grid") {
        s.plane = false;
        const auto &rows = j.at("z");
        if (!rows.is_array() || rows.size() != ny) invalid(label + ".z has the wrong row count");
        for (const auto &row : rows) {
            if (!row.is_array() || row.size() != nx) invalid(label + ".z has the wrong column count");
            for (const auto &v : row) s.z.push_back(number(v, label + ".z"));
        }
    } else invalid(label + ": unsupported surface type " + type);
    return s;
}
} // namespace

double Patch::value(XY p) const
{ const double u = p.x - x, v = p.y - y; return c + ax*u + ay*v + axy*u*v; }
XY Patch::gradient(XY p) const
{ return {ax + axy*(p.y-y), ay + axy*(p.x-x)}; }

LayerStack LayerStack::parse(const std::string &text)
{
    LayerStack out;
    try {
        const Json j = Json::parse(text);
        if (j.at("format") != "orca_nonplanar_layers") invalid("unknown format");
        if (!j.at("version").is_number_integer() || j.at("version") != 1) invalid("version must be 1");
        if (j.at("units") != "mm") invalid("units must be mm");
        if (j.at("coordinate_system") != "plate_world_xyz") invalid("coordinate_system must be plate_world_xyz");
        if (!j.at("plate_index").is_number_integer() || j.at("plate_index").get<int>() < 0)
            invalid("plate_index must be a nonnegative integer");
        out.plate_index = j.at("plate_index").get<int>();
        out.x = axis(j.at("grid").at("x"), "x");
        out.y = axis(j.at("grid").at("y"), "y");
        out.bed = surface(j.at("bed"), out.x.size(), out.y.size(), "bed");
        if (!out.bed.plane) invalid("bed must be a plane");
        const auto &l = j.at("limits");
        out.limits = {number(l.at("geometry_tolerance_mm"), "geometry_tolerance_mm"),
                      number(l.at("path_chord_tolerance_mm"), "path_chord_tolerance_mm"),
                      number(l.at("max_segment_length_mm"), "max_segment_length_mm"),
                      number(l.at("min_vertical_gap_mm"), "min_vertical_gap_mm"),
                      number(l.at("max_vertical_gap_mm"), "max_vertical_gap_mm"),
                      number(l.at("max_slope_degrees"), "max_slope_degrees")};
        const auto &lim = out.limits;
        if (lim.geometry_tolerance < 1e-5 || lim.chord_tolerance < 1e-5 || lim.max_segment_length <= 0 ||
            lim.min_gap <= 0 || lim.max_gap < lim.min_gap || lim.max_slope_degrees <= 0 || lim.max_slope_degrees >= 90)
            invalid("tolerances must be >= 0.00001 mm, lengths/gaps positive, min gap <= max gap, and slope in (0,90) degrees");
        const auto &layers = j.at("layers");
        if (!layers.is_array() || layers.empty()) invalid("layers must be a nonempty array");
        for (size_t k = 0; k < layers.size(); ++k) {
            const auto &entry = layers[k];
            if (!entry.at("id").is_number_integer() || entry.at("id").get<long long>() != static_cast<long long>(k))
                invalid("layer IDs must be consecutive from zero");
            out.layers.push_back(surface(entry.at("surface"), out.x.size(), out.y.size(), "layer " + std::to_string(k)));
        }
        const double slope = std::tan(lim.max_slope_degrees * std::acos(-1.) / 180.);
        for (size_t k = 0; k < out.layers.size(); ++k) {
            const Surface &s = out.layers[k], &prev = k ? out.layers[k-1] : out.bed;
            for (size_t j = 0; j+1 < out.y.size(); ++j)
                for (size_t i = 0; i+1 < out.x.size(); ++i) {
                    const auto patch = out.patch(s, i, j), below = out.patch(prev, i, j);
                    for (const XY p : {XY{out.x[i],out.y[j]}, XY{out.x[i+1],out.y[j]},
                                       XY{out.x[i],out.y[j+1]}, XY{out.x[i+1],out.y[j+1]}}) {
                        const double gap = patch.value(p) - below.value(p);
                        if (!std::isfinite(gap) || gap < lim.min_gap-1e-10 || gap > lim.max_gap+1e-10)
                            invalid("layer " + std::to_string(k) + ": vertical gap outside limits at grid cell " +
                                    std::to_string(i) + "," + std::to_string(j));
                        const XY g = patch.gradient(p);
                        // Norm of an affine gradient is convex; its maximum over a rectangle
                        // is attained at a corner, including both sides of grid seams.
                        if (!std::isfinite(g.x) || !std::isfinite(g.y) || std::hypot(g.x,g.y) > slope+1e-10)
                            invalid("layer " + std::to_string(k) + ": interpolated slope exceeds limit");
                    }
                }
        }
        out.source = j.dump();
    } catch (const Json::exception &e) { invalid(e.what()); }
    return out;
}

LayerStack LayerStack::read(const std::string &filename)
{
    std::ifstream file(std::filesystem::u8path(filename), std::ios::binary);
    if (!file) invalid("cannot open " + filename);
    return parse(std::string(std::istreambuf_iterator<char>(file), {}));
}

bool LayerStack::covers(XY p) const
{ return p.x >= x.front() && p.x <= x.back() && p.y >= y.front() && p.y <= y.back(); }

std::pair<size_t,size_t> LayerStack::cell(XY p) const
{
    if (!covers(p)) invalid("model or path outside shared grid at " + std::to_string(p.x) + "," + std::to_string(p.y));
    const size_t i = std::min(size_t(std::upper_bound(x.begin(),x.end(),p.x)-x.begin()-1),x.size()-2);
    const size_t j = std::min(size_t(std::upper_bound(y.begin(),y.end(),p.y)-y.begin()-1),y.size()-2);
    return {i,j};
}

Patch LayerStack::patch(const Surface &s, size_t i, size_t j) const
{
    if (s.plane) return {x[i],y[j],s.a*x[i]+s.b*y[j]+s.c,s.a,s.b,0.};
    const size_t n = x.size();
    const double dx = x[i+1]-x[i], dy = y[j+1]-y[j];
    const double z00 = s.z[j*n+i], z10 = s.z[j*n+i+1], z01 = s.z[(j+1)*n+i], z11 = s.z[(j+1)*n+i+1];
    return {x[i],y[j],z00,(z10-z00)/dx,(z01-z00)/dy,(z11-z10-z01+z00)/(dx*dy)};
}

double LayerStack::value(const Surface &s, XY p) const
{ const auto [i,j] = cell(p); return patch(s,i,j).value(p); }
XY LayerStack::gradient(const Surface &s, XY p) const
{ const auto [i,j] = cell(p); return patch(s,i,j).gradient(p); }
std::pair<double,double> LayerStack::bounds(const Surface &s) const
{
    if (!s.plane) {
        const auto [lo,hi] = std::minmax_element(s.z.begin(),s.z.end());
        return {*lo,*hi};
    }
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (double py : {y.front(),y.back()}) for (double px : {x.front(),x.back()}) {
        const double z = value(s,{px,py}); lo = std::min(lo,z); hi = std::max(hi,z);
    }
    return {lo,hi};
}
double LayerStack::vertical_gap(size_t k, XY p) const
{ return value(layers[k],p)-value(k ? layers[k-1] : bed,p); }
double LayerStack::normal_spacing(size_t k, XY p) const
{
    const Surface &previous = k ? layers[k-1] : bed;
    const XY g = gradient(layers[k],p);
    const double norm = std::sqrt(1.+g.x*g.x+g.y*g.y);
    const XY direction {g.x/norm,g.y/norm}; // along the downward top normal
    const double dz = -1./norm, top = value(layers[k],p);
    if (g.x == 0. && g.y == 0.) return vertical_gap(k,p);
    double end = std::numeric_limits<double>::infinity();
    for (const auto [coordinate,d,lo,hi] : {std::array<double,4>{p.x,direction.x,x.front(),x.back()},
                                           std::array<double,4>{p.y,direction.y,y.front(),y.back()}})
        if (d != 0.) end = std::min(end,((d > 0. ? hi : lo)-coordinate)/d);
    // The ray can never meet a preceding positive-height sheet after passing
    // below its global lower bound. This also bounds nearly vertical rays.
    end = std::min(end,std::max(0.,(top-bounds(previous).first)/(-dz))+limits.max_gap);
    std::vector<double> cuts {0.,end};
    if (direction.x != 0.) for (double px : x) {
        const double t = (px-p.x)/direction.x; if (t > 0 && t < end) cuts.push_back(t);
    }
    if (direction.y != 0.) for (double py : y) {
        const double t = (py-p.y)/direction.y; if (t > 0 && t < end) cuts.push_back(t);
    }
    std::sort(cuts.begin(),cuts.end());
    for (size_t n = 1; n < cuts.size(); ++n) {
        const double t0 = cuts[n-1], span = cuts[n]-t0;
        if (span <= 0.) continue;
        const XY mid {p.x+direction.x*(t0+.5*span),p.y+direction.y*(t0+.5*span)};
        const auto [i,j] = cell(mid);
        const Patch patch = this->patch(previous,i,j);
        auto f = [&](double u) {
            const double t = t0+u*span;
            return top+dz*t-patch.value({p.x+direction.x*t,p.y+direction.y*t});
        };
        const double c = f(0.), f1 = f(1.), fm = f(.5);
        const double a = 2.*(f1+c-2.*fm), b = f1-c-a;
        std::vector<double> roots;
        if (std::abs(a) < 1e-12) { if (b != 0.) roots.push_back(-c/b); }
        else {
            const double discriminant = b*b-4.*a*c;
            if (discriminant >= 0.) {
                const double root = std::sqrt(discriminant);
                // Stable quadratic evaluation, including very unequal roots.
                const double q = -.5*(b+std::copysign(root,b));
                roots.push_back(q/a);
                if (q != 0.) roots.push_back(c/q);
            }
        }
        std::sort(roots.begin(),roots.end());
        for (double u : roots) if (u >= -1e-10 && u <= 1.+1e-10 && t0+u*span > 0.) return t0+u*span;
    }
    invalid("downward normal does not reach the preceding surface within the shared grid; enlarge the domain or redesign the stack");
}

} // namespace Slic3r::Nonplanar
