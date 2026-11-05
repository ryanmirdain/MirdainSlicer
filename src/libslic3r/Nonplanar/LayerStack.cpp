#include "LayerStack.hpp"
#include "../I18N.hpp"
#include "../format.hpp"

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
{ throw std::invalid_argument(Slic3r::format(_u8L("Nonplanar layers: %1%"),message)); }

double number(const Json &j, const std::string &label)
{
    if (!j.is_number() || !std::isfinite(j.get<double>())) invalid(Slic3r::format(_u8L("%1% must be finite"),label));
    return j.get<double>();
}

std::vector<double> axis(const Json &j, const std::string &name)
{
    if (!j.is_array() || j.size() < 2) invalid(Slic3r::format(_u8L("grid.%1% needs at least two values"),name));
    std::vector<double> out;
    for (const auto &v : j) {
        const double n = number(v, "grid." + name);
        // Clipper coordinates use 1e-6 mm. Larger coordinates are not printer coordinates.
        if (std::abs(n) > 1e7) invalid(_u8L("grid coordinates exceed the supported physical range"));
        if (!out.empty() && n <= out.back()) invalid(Slic3r::format(_u8L("grid.%1% must be strictly increasing"),name));
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
        if (!rows.is_array() || rows.size() != ny) invalid(Slic3r::format(_u8L("%1%.z has the wrong row count"),label));
        for (const auto &row : rows) {
            if (!row.is_array() || row.size() != nx) invalid(Slic3r::format(_u8L("%1%.z has the wrong column count"),label));
            for (const auto &v : row) s.z.push_back(number(v, label + ".z"));
        }
        const auto bounds = std::minmax_element(s.z.begin(),s.z.end());
        s.min_z = *bounds.first; s.max_z = *bounds.second;
    } else invalid(Slic3r::format(_u8L("%1%: unsupported surface type %2%"),label,type));
    return s;
}
// Convert a clamped cubic B-spline to Bezier spans using knot insertion.
// Controls at Greville sites are convex combinations of the supplied samples.
// Knot insertion is also convex, so every span uses positive common weights.
std::vector<std::vector<double>> spline_weights(const std::vector<double> &axis)
{
    std::vector<double> knots(4,axis.front());
    knots.insert(knots.end(),axis.begin()+1,axis.end()-1);
    knots.insert(knots.end(),4,axis.back());
    std::vector<std::vector<double>> controls(knots.size()-4,std::vector<double>(axis.size(),0.));
    for (size_t n = 0; n < controls.size(); ++n) {
        const double p = std::clamp((knots[n+1]+knots[n+2]+knots[n+3])/3.,axis.front(),axis.back());
        const size_t i = std::min(size_t(std::upper_bound(axis.begin(),axis.end(),p)-axis.begin()-1),axis.size()-2);
        const double t = std::clamp((p-axis[i])/(axis[i+1]-axis[i]),0.,1.);
        controls[n][i] = 1.-t; controls[n][i+1] = t;
    }
    for (size_t n = 1; n+1 < axis.size(); ++n) for (size_t repeat = 0; repeat < 2; ++repeat) {
        const double u = axis[n];
        const size_t k = std::upper_bound(knots.begin(),knots.end(),u)-knots.begin()-1;
        const size_t multiplicity = std::count(knots.begin(),knots.end(),u);
        std::vector<std::vector<double>> next(controls.size()+1,std::vector<double>(axis.size()));
        for (size_t i = 0; i <= k-3; ++i) next[i] = controls[i];
        for (size_t i = k-multiplicity; i < controls.size(); ++i) next[i+1] = controls[i];
        for (size_t i = k-2; i <= k-multiplicity; ++i) {
            const double alpha = (u-knots[i])/(knots[i+3]-knots[i]);
            for (size_t j = 0; j < axis.size(); ++j)
                next[i][j] = alpha*controls[i][j]+(1.-alpha)*controls[i-1][j];
        }
        controls = std::move(next);
        knots.insert(knots.begin()+k+1,u);
    }
    return controls;
}
std::array<double,4> basis(double t)
{ const double s = 1.-t; return {s*s*s,3.*s*s*t,3.*s*t*t,t*t*t}; }
std::array<double,4> derivative(double t)
{ const double s = 1.-t; return {-3.*s*s,3.*s*s-6.*s*t,6.*s*t-3.*t*t,3.*t*t}; }
std::array<double,4> restrict_curve(const std::array<double,4> &control, double a, double b)
{
    std::array<double,4> out;
    for (size_t n = 0; n < 4; ++n) {
        auto q = control;
        for (size_t step = 0; step < 3; ++step) {
            const double t = step < 3-n ? a : b;
            for (size_t i = 0; i < 3-step; ++i) q[i] = (1.-t)*q[i]+t*q[i+1];
        }
        out[n] = q[0];
    }
    return out;
}
} // namespace

double Patch::value(XY p) const
{
    const auto bx = basis((p.x-x)/dx), by = basis((p.y-y)/dy);
    double value = 0.;
    for (size_t j = 0; j < 4; ++j) for (size_t i = 0; i < 4; ++i)
        value += control[j*4+i]*bx[i]*by[j];
    return value;
}
XY Patch::gradient(XY p) const
{
    const double u = (p.x-x)/dx, v = (p.y-y)/dy;
    const auto bx = basis(u), by = basis(v), gx = derivative(u), gy = derivative(v);
    XY g {0.,0.};
    for (size_t j = 0; j < 4; ++j) for (size_t i = 0; i < 4; ++i) {
        g.x += control[j*4+i]*gx[i]*by[j]/dx;
        g.y += control[j*4+i]*bx[i]*gy[j]/dy;
    }
    return g;
}
std::pair<double,double> Patch::bounds() const
{
    const auto bounds = std::minmax_element(control.begin(),control.end());
    return {*bounds.first,*bounds.second};
}
double Patch::linear_error() const
{
    double xx = 0., yy = 0., xy = 0.;
    for (size_t j = 0; j < 4; ++j) for (size_t i = 0; i < 2; ++i)
        xx = std::max(xx,6.*std::abs(control[j*4+i+2]-2.*control[j*4+i+1]+control[j*4+i]));
    for (size_t j = 0; j < 2; ++j) for (size_t i = 0; i < 4; ++i)
        yy = std::max(yy,6.*std::abs(control[(j+2)*4+i]-2.*control[(j+1)*4+i]+control[j*4+i]));
    for (size_t j = 0; j < 3; ++j) for (size_t i = 0; i < 3; ++i)
        xy = std::max(xy,9.*std::abs(control[(j+1)*4+i+1]-control[(j+1)*4+i]-control[j*4+i+1]+control[j*4+i]));
    // Taylor's remainder on any triangle in the normalized unit cell.
    return .5*(xx+2.*xy+yy);
}
std::array<double,7> Patch::line(XY a, XY b) const
{
    std::array<double,16> q;
    for (size_t j = 0; j < 4; ++j) {
        const auto row = restrict_curve({control[j*4],control[j*4+1],control[j*4+2],control[j*4+3]},(a.x-x)/dx,(b.x-x)/dx);
        for (size_t i = 0; i < 4; ++i) q[j*4+i] = row[i];
    }
    for (size_t i = 0; i < 4; ++i) {
        const auto column = restrict_curve({q[i],q[4+i],q[8+i],q[12+i]},(a.y-y)/dy,(b.y-y)/dy);
        for (size_t j = 0; j < 4; ++j) q[j*4+i] = column[j];
    }
    constexpr std::array<double,4> c3 {1.,3.,3.,1.};
    constexpr std::array<double,7> c6 {1.,6.,15.,20.,15.,6.,1.};
    std::array<double,7> curve {};
    for (size_t j = 0; j < 4; ++j) for (size_t i = 0; i < 4; ++i)
        curve[i+j] += q[j*4+i]*c3[i]*c3[j]/c6[i+j];
    return curve;
}

LayerStack LayerStack::parse(const std::string &text)
{
    LayerStack out;
    try {
        Json j = Json::parse(text);
        if (j.at("format") != "orca_nonplanar_layers") invalid(_u8L("unknown format"));
        if (!j.at("version").is_number_integer() || j.at("version") != 1) invalid(_u8L("version must be 1"));
        if (j.at("units") != "mm") invalid(_u8L("units must be mm"));
        if (j.at("coordinate_system") != "plate_world_xyz") invalid(_u8L("coordinate_system must be plate_world_xyz"));
        if (!j.at("plate_index").is_number_integer() || j.at("plate_index") < 0 ||
            j.at("plate_index") > std::numeric_limits<int>::max())
            invalid(_u8L("plate_index must be a nonnegative integer within the supported plate-number range"));
        out.plate_index = j.at("plate_index").get<int>();
        out.x = axis(j.at("grid").at("x"), "x");
        out.y = axis(j.at("grid").at("y"), "y");
        out.bed = surface(j.at("bed"), out.x.size(), out.y.size(), "bed");
        const auto &l = j.at("limits");
        out.limits = {number(l.at("geometry_tolerance_mm"), "geometry_tolerance_mm"),
                      number(l.at("path_chord_tolerance_mm"), "path_chord_tolerance_mm"),
                      number(l.at("max_segment_length_mm"), "max_segment_length_mm"),
                      number(l.at("min_vertical_gap_mm"), "min_vertical_gap_mm"),
                      number(l.at("max_vertical_gap_mm"), "max_vertical_gap_mm"),
                      number(l.at("max_slope_degrees"), "max_slope_degrees")};
        const auto &lim = out.limits;
        if (lim.geometry_tolerance <= 0. || lim.chord_tolerance <= 0. || lim.max_segment_length <= 0.)
            invalid(_u8L("geometry/chord tolerances and maximum segment length must be positive"));
        const auto &layers = j.at("layers");
        if (!layers.is_array() || layers.empty()) invalid(_u8L("layers must be a nonempty array"));
        for (size_t k = 0; k < layers.size(); ++k) {
            const auto &entry = layers[k];
            if (!entry.at("id").is_number_integer() || entry.at("id").get<long long>() != static_cast<long long>(k))
                invalid(_u8L("layer IDs must be consecutive from zero"));
            out.layers.push_back(surface(entry.at("surface"), out.x.size(), out.y.size(), "layer " + std::to_string(k)));
        }
        auto sample = [&](const Surface &s, size_t i, size_t j) {
            return s.plane ? s.a*out.x[i]+s.b*out.y[j]+s.c : s.z[j*out.x.size()+i];
        };
        for (size_t k = 0; k < out.layers.size(); ++k)
            for (size_t j = 0; j < out.y.size(); ++j) for (size_t i = 0; i < out.x.size(); ++i) {
                const double gap = sample(out.layers[k],i,j)-sample(k ? out.layers[k-1] : out.bed,i,j);
                if (!std::isfinite(gap) || gap <= 0.)
                    invalid(Slic3r::format(_u8L("Layer %1%: surfaces must be strictly ordered at grid point %2%,%3%"),k,i,j));
            }
        const auto wx = spline_weights(out.x), wy = spline_weights(out.y);
        auto build = [&](Surface &s) {
            std::vector<double> rows(out.y.size()*wx.size(),0.);
            for (size_t v = 0; v < out.y.size(); ++v) for (size_t a = 0; a < wx.size(); ++a)
                for (size_t u = 0; u < out.x.size(); ++u) if (wx[a][u] != 0.)
                    rows[v*wx.size()+a] += wx[a][u]*sample(s,u,v);
            std::vector<double> controls(wx.size()*wy.size(),0.);
            for (size_t b = 0; b < wy.size(); ++b) for (size_t a = 0; a < wx.size(); ++a)
                for (size_t v = 0; v < out.y.size(); ++v) if (wy[b][v] != 0.)
                    controls[b*wx.size()+a] += wy[b][v]*rows[v*wx.size()+a];
            const auto bounds = std::minmax_element(controls.begin(),controls.end());
            s.min_z = *bounds.first; s.max_z = *bounds.second;
            for (size_t j = 0; j+1 < out.y.size(); ++j) for (size_t i = 0; i+1 < out.x.size(); ++i) {
                Patch patch {out.x[i],out.y[j],out.x[i+1]-out.x[i],out.y[j+1]-out.y[j]};
                for (size_t b = 0; b < 4; ++b) for (size_t a = 0; a < 4; ++a)
                    patch.control[b*4+a] = controls[(j*3+b)*wx.size()+i*3+a];
                s.patches.push_back(std::move(patch));
            }
        };
        build(out.bed);
        for (auto &s : out.layers) build(s);
        j["surface_interpolation"] = "ordered_cubic_bspline_v1";
        out.source = j.dump();
    } catch (const Json::exception &e) { invalid(e.what()); }
    return out;
}

LayerStack LayerStack::read(const std::string &filename)
{
    std::ifstream file(std::filesystem::u8path(filename), std::ios::binary);
    if (!file) invalid(Slic3r::format(_u8L("Cannot open %1%"),filename));
    return parse(std::string(std::istreambuf_iterator<char>(file), {}));
}

bool LayerStack::covers(XY p) const
{ return p.x >= x.front() && p.x <= x.back() && p.y >= y.front() && p.y <= y.back(); }

std::pair<size_t,size_t> LayerStack::cell(XY p) const
{
    if (!covers(p)) invalid(Slic3r::format(_u8L("Model or path outside shared grid at %1%,%2%"),p.x,p.y));
    const size_t i = std::min(size_t(std::upper_bound(x.begin(),x.end(),p.x)-x.begin()-1),x.size()-2);
    const size_t j = std::min(size_t(std::upper_bound(y.begin(),y.end(),p.y)-y.begin()-1),y.size()-2);
    return {i,j};
}

Patch LayerStack::patch(const Surface &s, size_t i, size_t j) const
{ return s.patches[j*(x.size()-1)+i]; }

double LayerStack::value(const Surface &s, XY p) const
{ const auto [i,j] = cell(p); return patch(s,i,j).value(p); }
XY LayerStack::gradient(const Surface &s, XY p) const
{ const auto [i,j] = cell(p); return patch(s,i,j).gradient(p); }
std::pair<double,double> LayerStack::bounds(const Surface &s) const
{ return {s.min_z,s.max_z}; }
double LayerStack::vertical_gap(size_t k, XY p) const
{ return value(layers[k],p)-value(k ? layers[k-1] : bed,p); }
double LayerStack::normal_spacing(size_t k, XY p) const
{
    const XY g = gradient(layers[k],p);
    return vertical_gap(k,p)/std::sqrt(1.+g.x*g.x+g.y*g.y);
}

} // namespace Slic3r::Nonplanar
