#include "Slicer.hpp"
#include <clipper2/clipper.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Slic3r::Nonplanar {
namespace {
using namespace Clipper2Lib;
constexpr double scale = 1e6;
constexpr double boundary_epsilon = 1e-8;
[[noreturn]] void invalid(const std::string &s) { throw std::invalid_argument("Nonplanar slicing: " + s); }
XY xy(XYZ p) { return {p.x,p.y}; }
XYZ mix(XYZ a, XYZ b, double t)
{ return {a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t,a.z+(b.z-a.z)*t}; }
double distance(XYZ a, XYZ b)
{ return std::hypot(std::hypot(a.x-b.x,a.y-b.y),a.z-b.z); }
double area(XYZ a, XYZ b, XYZ c)
{ return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x); }
Point64 point(XY p) { return {std::llround(p.x*scale),std::llround(p.y*scale)}; }
XY real(Point64 p) { return {p.x/scale,p.y/scale}; }
Paths64 unite(const Paths64 &p) { return Union(p,FillRule::NonZero); }
Paths64 intersect(const Paths64 &a, const Paths64 &b) { return Intersect(a,b,FillRule::NonZero); }
Paths64 difference(const Paths64 &a, const Paths64 &b) { return Difference(a,b,FillRule::NonZero); }
Paths64 inset(const Paths64 &p, double d, double tolerance)
{ return InflatePaths(p,-d*scale,Clipper2Lib::JoinType::Round,Clipper2Lib::EndType::Polygon,2.,tolerance*scale/8.); }

// Clip in original XYZ against XY cell edges. No model deformation or exported mesh.
std::vector<XYZ> clip_axis(const std::vector<XYZ> &polygon, int axis, double value, bool lower)
{
    std::vector<XYZ> out;
    if (polygon.empty()) return out;
    auto f = [=](XYZ p) { return (axis == 0 ? p.x : p.y)-value; };
    XYZ a = polygon.back();
    double fa = f(a);
    for (XYZ b : polygon) {
        const double fb = f(b);
        const bool ina = lower ? fa >= 0 : fa <= 0, inb = lower ? fb >= 0 : fb <= 0;
        if (ina != inb) out.push_back(mix(a,b,fa/(fa-fb)));
        if (inb) out.push_back(b);
        a = b; fa = fb;
    }
    return out;
}

class Intersection {
public:
    Intersection(const LayerStack &s, const std::function<void()> &c) : stack(s), cancel(c) {}

    void require_coverage(const Object &object)
    {
        Paths64 uncovered;
        for (const auto &t : object.triangles)
            face(t,stack.layers.back(),stack.limits.geometry_tolerance,uncovered,true);
        if (!uncovered.empty()) invalid("final supplied surface lies below the original mesh between vertices");
    }

    // Positive/negative projected face winding implements an upward vertical ray
    // through the original closed oriented solid. Vertical faces have zero area.
    Paths64 material(const Object &object, const Surface &surface, double shift)
    {
        Paths64 polygons;
        for (const auto &t : object.triangles)
            face(t,surface,shift,polygons,false);
        return Union(polygons,FillRule::Positive);
    }

    Paths64 band(const Object &object, size_t k)
    {
        const Surface &top = stack.layers[k], &bottom = k ? stack.layers[k-1] : stack.bed;
        Paths64 region = intersect(material(object,top,-boundary_epsilon), material(object,bottom,boundary_epsilon));
        // Endpoint containment alone misses cavities lying wholly between the
        // sheets. Remove every projected boundary strictly inside the band.
        Paths64 boundary;
        for (const auto &t : object.triangles) {
            Paths64 above_bottom, above_top;
            face(t,bottom,boundary_epsilon,above_bottom,true);
            face(t,top,-boundary_epsilon,above_top,true);
            const auto between = difference(unite(above_bottom),unite(above_top));
            boundary.insert(boundary.end(),between.begin(),between.end());
        }
        return difference(region,unite(boundary));
    }

private:
    const LayerStack &stack;
    const std::function<void()> &cancel;
    size_t nodes = 0;

    void append(const std::vector<XYZ> &vertices, Paths64 &out, bool positive)
    {
        Path64 polygon;
        for (XYZ p : vertices) polygon.push_back(point(xy(p)));
        if (polygon.size() < 3 || std::abs(Area(polygon)) < 1.) return;
        if (positive && !IsPositive(polygon)) std::reverse(polygon.begin(),polygon.end());
        out.push_back(std::move(polygon));
    }

    void refine(std::array<XYZ,3> t, const Patch &h, double shift, Paths64 &out, bool positive, unsigned depth = 0)
    {
        if ((++nodes & 4095) == 0) cancel();
        if (nodes > 20000000) invalid("intersection refinement budget exceeded; simplify the layer grid or relax geometry_tolerance_mm");
        auto f = [&](XYZ p) { return p.z-h.value(xy(p))-shift; };
        std::array<double,3> v {f(t[0]),f(t[1]),f(t[2])};
        std::array<XYZ,3> mid {mix(t[0],t[1],.5),mix(t[1],t[2],.5),mix(t[2],t[0],.5)};
        std::array<double,3> edge {};
        double lo = *std::min_element(v.begin(),v.end()), hi = *std::max_element(v.begin(),v.end());
        double error = 0.;
        for (size_t e = 0; e < 3; ++e) {
            // Degree-two Bernstein control values bound the full quadratic on
            // this triangle, including same-sign vertex crossings and tangencies.
            const double linear = .5*(v[e]+v[(e+1)%3]);
            edge[e] = 2.*f(mid[e])-linear;
            lo = std::min(lo,edge[e]); hi = std::max(hi,edge[e]);
            error = std::max(error,std::abs(edge[e]-linear));
        }
        if (lo >= 0.) { append({t.begin(),t.end()},out,positive); return; }
        if (hi < 0.) return;
        size_t longest = 0;
        double length = 0.;
        for (size_t e = 0; e < 3; ++e) {
            const double d = std::hypot(t[e].x-t[(e+1)%3].x,t[e].y-t[(e+1)%3].y);
            if (d > length) { length = d; longest = e; }
        }
        if (error < 1e-12 || length <= stack.limits.geometry_tolerance*.5) {
            // Linear clipping on a bounded local patch. Positive and negative
            // winding use opposite conservative sides of the residual bound.
            const double margin = positive || area(t[0],t[1],t[2]) > 0. ? error : -error;
            std::vector<XYZ> polygon;
            XYZ a = t.back(); double fa = v.back()-margin;
            for (size_t e = 0; e < 3; ++e) {
                XYZ b = t[e]; const double fb = v[e]-margin;
                if ((fa >= 0) != (fb >= 0)) polygon.push_back(mix(a,b,fa/(fa-fb)));
                if (fb >= 0) polygon.push_back(b);
                a = b; fa = fb;
            }
            append(polygon,out,positive);
            return;
        }
        if (depth == 48) invalid("intersection did not converge within geometry tolerance");
        const size_t a = longest, b = (a+1)%3, c = (a+2)%3;
        refine({t[a],mid[a],t[c]},h,shift,out,positive,depth+1);
        refine({mid[a],t[b],t[c]},h,shift,out,positive,depth+1);
    }

    void face(const Triangle &t, const Surface &surface, double shift, Paths64 &out, bool positive)
    {
        if (std::abs(area(t.vertices[0],t.vertices[1],t.vertices[2])) < 1e-16) return;
        double xmin = t.vertices[0].x, xmax = xmin, ymin = t.vertices[0].y, ymax = ymin;
        for (XYZ p : t.vertices) {
            xmin = std::min(xmin,p.x); xmax = std::max(xmax,p.x);
            ymin = std::min(ymin,p.y); ymax = std::max(ymax,p.y);
        }
        const auto [i0,j0] = stack.cell({xmin,ymin});
        const auto [i1,j1] = stack.cell({xmax,ymax});
        for (size_t j = j0; j <= j1; ++j) for (size_t i = i0; i <= i1; ++i) {
            std::vector<XYZ> p(t.vertices.begin(),t.vertices.end());
            p = clip_axis(p,0,stack.x[i],true); p = clip_axis(p,0,stack.x[i+1],false);
            p = clip_axis(p,1,stack.y[j],true); p = clip_axis(p,1,stack.y[j+1],false);
            for (size_t n = 1; n+1 < p.size(); ++n)
                refine({p[0],p[n],p[n+1]},stack.patch(surface,i,j),shift,out,positive);
        }
    }
};

class Planner {
public:
    Planner(const LayerStack &s, const Settings &p) : stack(s), settings(p) {}

    Path path(const Path64 &points, bool closed, size_t k, Role role)
    {
        Path out;
        out.role = role;
        const size_t seam = closed ? size_t(std::min_element(points.begin(),points.end(),[](Point64 a, Point64 b) {
            return a.x == b.x ? a.y < b.y : a.x < b.x;
        })-points.begin()) : 0;
        const XY first = real(points[seam]);
        out.start = {first.x,first.y,stack.value(stack.layers[k],first)};
        for (size_t n = 1; n < points.size()+size_t(closed); ++n) {
            const XY a = real(points[(seam+n-1)%points.size()]), b = real(points[(seam+n)%points.size()]);
            // Split at every grid seam before curvature/flow subdivision. A
            // bilinear segment is quadratic within a cell, never across cells.
            std::vector<double> cuts {0.,1.};
            if (a.x != b.x) for (double x : stack.x) { const double t = (x-a.x)/(b.x-a.x); if (t > 0 && t < 1) cuts.push_back(t); }
            if (a.y != b.y) for (double y : stack.y) { const double t = (y-a.y)/(b.y-a.y); if (t > 0 && t < 1) cuts.push_back(t); }
            std::sort(cuts.begin(),cuts.end());
            cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            for (size_t c = 1; c < cuts.size(); ++c) {
                auto at = [&](double t) { return XY{a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t}; };
                segment(out,at(cuts[c-1]),at(cuts[c]),k,0);
            }
        }
        return out;
    }

private:
    const LayerStack &stack;
    const Settings &settings;

    // Integrate the actual bilinear band over the complete projected rectangular
    // swept cross-section, splitting the strip at cell boundaries. Polygon
    // moments integrate 1, x, y and xy exactly; this includes transverse changes
    // in thickness which a centerline-only length correction would miss.
    double strip_volume(XY a, XY b, size_t k, double &preview_error) const
    {
        const double ds = std::hypot(b.x-a.x,b.y-a.y), half = settings.width*.5;
        const XY normal {-(b.y-a.y)/ds*half,(b.x-a.x)/ds*half};
        const std::vector<XYZ> strip {{a.x+normal.x,a.y+normal.y,0.},{b.x+normal.x,b.y+normal.y,0.},
                                     {b.x-normal.x,b.y-normal.y,0.},{a.x-normal.x,a.y-normal.y,0.}};
        double xmin = a.x, xmax = a.x, ymin = a.y, ymax = a.y;
        for (XYZ p : strip) { xmin = std::min(xmin,p.x); xmax = std::max(xmax,p.x); ymin = std::min(ymin,p.y); ymax = std::max(ymax,p.y); }
        const auto [i0,j0] = stack.cell({xmin,ymin});
        const auto [i1,j1] = stack.cell({xmax,ymax});
        double volume = 0.;
        const Surface &below = k ? stack.layers[k-1] : stack.bed;
        const XY mid {(a.x+b.x)*.5,(a.y+b.y)*.5};
        const XY grad = stack.gradient(stack.layers[k],mid);
        const double za = stack.value(stack.layers[k],a), zb = stack.value(stack.layers[k],b);
        const double gap = stack.vertical_gap(k,mid);
        const XY tangent {(b.x-a.x)/ds,(b.y-a.y)/ds}, right {tangent.y,-tangent.x};
        const double cross_slope = grad.x*right.x+grad.y*right.y;
        const XY slope {tangent.x*(zb-za)/ds+right.x*cross_slope,
                        tangent.y*(zb-za)/ds+right.y*cross_slope};
        preview_error = 0.;
        for (size_t j = j0; j <= j1; ++j) for (size_t i = i0; i <= i1; ++i) {
            auto polygon = clip_axis(strip,0,stack.x[i],true); polygon = clip_axis(polygon,0,stack.x[i+1],false);
            polygon = clip_axis(polygon,1,stack.y[j],true); polygon = clip_axis(polygon,1,stack.y[j+1],false);
            if (polygon.size() < 3) continue;
            const Patch top = stack.patch(stack.layers[k],i,j), bottom = stack.patch(below,i,j);
            // Bound the linear preview's complete top AND bottom, including
            // transverse curvature. A bead that cannot be represented within
            // tolerance is rejected rather than displayed as a planar bead.
            for (int side = 0; side < 2; ++side) {
                const Patch &sheet = side ? bottom : top;
                auto residual = [&](XYZ p) {
                    return sheet.value(xy(p))-(za+slope.x*(p.x-a.x)+slope.y*(p.y-a.y)-(side ? gap : 0.));
                };
                for (size_t n = 1; n+1 < polygon.size(); ++n) {
                    const std::array<XYZ,3> tri {polygon[0],polygon[n],polygon[n+1]};
                    for (size_t v = 0; v < 3; ++v) {
                        const double f0 = residual(tri[v]), f1 = residual(tri[(v+1)%3]);
                        const double control = 2.*residual(mix(tri[v],tri[(v+1)%3],.5))-.5*(f0+f1);
                        preview_error = std::max({preview_error,std::abs(f0),std::abs(control)});
                    }
                }
            }
            double area = 0., mx = 0., my = 0., mxy = 0.;
            for (size_t n = 0; n < polygon.size(); ++n) {
                const XY p {polygon[n].x-top.x,polygon[n].y-top.y};
                const XYZ next = polygon[(n+1)%polygon.size()];
                const XY q {next.x-top.x,next.y-top.y};
                const double cross = p.x*q.y-q.x*p.y;
                area += cross*.5; mx += (p.x+q.x)*cross/6.; my += (p.y+q.y)*cross/6.;
                mxy += (2.*p.x*p.y+p.x*q.y+q.x*p.y+2.*q.x*q.y)*cross/24.;
            }
            volume += std::abs((top.c-bottom.c)*area+(top.ax-bottom.ax)*mx+(top.ay-bottom.ay)*my+(top.axy-bottom.axy)*mxy);
        }
        return volume;
    }

    void segment(Path &out, XY a, XY b, size_t k, unsigned depth)
    {
        const XY mid {(a.x+b.x)*.5,(a.y+b.y)*.5};
        const XYZ p {a.x,a.y,stack.value(stack.layers[k],a)}, q {b.x,b.y,stack.value(stack.layers[k],b)};
        const double zm = stack.value(stack.layers[k],mid), gap = stack.vertical_gap(k,mid);
        const double ga = stack.vertical_gap(k,a), gb = stack.vertical_gap(k,b);
        if (std::max({ga,gb,gap}) > settings.nozzle_diameter)
            invalid("local vertical gap exceeds the nozzle diameter");
        const double chord_error = std::abs(zm-(p.z+q.z)*.5);
        // Include the cross-width surface variation in the modeled bead top:
        // the emitted centerline is on H; the ideal band envelope follows H.
        if (distance(p,q) > stack.limits.max_segment_length-.00018 || chord_error > stack.limits.chord_tolerance*.5 ||
            std::max(std::abs(ga-gap),std::abs(gb-gap)) > .01*gap) {
            if (depth == 40) invalid("path/flow segmentation did not converge");
            segment(out,a,mid,k,depth+1); segment(out,mid,b,k,depth+1); return;
        }
        const double ds = std::hypot(b.x-a.x,b.y-a.y);
        if (ds < 1e-10) return;
        const XY grad = stack.gradient(stack.layers[k],mid);
        const double cross_slope = (-grad.x*(b.y-a.y)+grad.y*(b.x-a.x))/ds;
        double preview_error;
        const double volume = strip_volume(a,b,k,preview_error);
        if (preview_error > stack.limits.geometry_tolerance*.5) {
            if (depth >= 24) invalid("line width cannot represent transverse surface curvature/thickness within geometry tolerance; use a narrower line or a larger tolerance");
            segment(out,a,mid,k,depth+1); segment(out,mid,b,k,depth+1); return;
        }
        const double normal = stack.normal_spacing(k,mid);
        if ((settings.min_normal_height > 0. && normal < settings.min_normal_height) ||
            (settings.max_normal_height > 0. && normal > settings.max_normal_height))
            invalid("local normal spacing exceeds the nozzle profile's layer-height limits");
        out.segments.push_back({q,volume,gap,normal,settings.width*std::sqrt(1.+cross_slope*cross_slope),-cross_slope});
    }
};

Paths64 infill(const Paths64 &region, const LayerStack &stack, double spacing, double angle)
{
    Paths64 lines;
    if (region.empty()) return lines;
    const double radians = angle*std::acos(-1.)/180.;
    const XY tangent {std::cos(radians),std::sin(radians)}, normal {-tangent.y,tangent.x};
    double low = std::numeric_limits<double>::infinity(), high = -low, start = low, end = -low;
    for (double x : {stack.x.front(),stack.x.back()}) for (double y : {stack.y.front(),stack.y.back()}) {
        const double n = x*normal.x+y*normal.y, t = x*tangent.x+y*tangent.y;
        low = std::min(low,n); high = std::max(high,n); start = std::min(start,t); end = std::max(end,t);
    }
    for (double d = std::ceil(low/spacing)*spacing; d <= high; d += spacing) {
        lines.push_back({point({tangent.x*start+normal.x*d,tangent.y*start+normal.y*d}),
                         point({tangent.x*end+normal.x*d,tangent.y*end+normal.y*d})});
    }
    Clipper64 clipper;
    clipper.AddOpenSubject(lines); clipper.AddClip(region);
    Paths64 closed, open;
    clipper.Execute(Clipper2Lib::ClipType::Intersection,FillRule::NonZero,closed,open);
    return open;
}
} // namespace

std::vector<Layer> slice(const LayerStack &stack, const std::vector<Object> &objects, const Settings &settings,
                         const std::function<void()> &cancel)
{
    if (settings.walls < 1 || settings.width <= 0 || settings.density < 0 || settings.density > 1)
        invalid("at least one wall, a positive line width and density in [0,1] are required");
    const double tol = stack.limits.geometry_tolerance;
    for (const Object &object : objects) for (const Triangle &t : object.triangles) for (XYZ p : t.vertices) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) invalid("model contains nonfinite coordinates");
        if (!stack.covers(xy(p))) invalid("shared grid does not cover the placed model footprint");
        if (p.z < stack.value(stack.bed,xy(p))-tol) invalid("model extends below the supplied bed");
        if (p.z > settings.max_z+tol) invalid("model exceeds machine Z height");
        if (p.z > stack.value(stack.layers.back(),xy(p))+tol) invalid("layer stack does not cover the complete model height");
    }
    std::vector<Layer> output(stack.layers.size());
    std::vector<Paths64> occupied(stack.layers.size());
    Planner planner(stack,settings);
    for (const Object &object : objects) {
        Intersection(stack,cancel).require_coverage(object);
        std::vector<Paths64> bands;
        for (size_t k = 0; k < stack.layers.size(); ++k) {
            cancel();
            Intersection engine(stack,cancel);
            bands.push_back(engine.band(object,k));
            if (!intersect(occupied[k],bands.back()).empty()) invalid("overlapping printable objects are unsupported");
            occupied[k].insert(occupied[k].end(),bands.back().begin(),bands.back().end());
        }
        for (size_t k = 0; k < bands.size(); ++k) {
            cancel();
            Layer &layer = output[k];
            const double margin = settings.width*.5+tol;
            Paths64 safe = inset(bands[k],margin,tol);
            if (k > 0 && !difference(safe,bands[k-1]).empty())
                invalid("unsupported overhang or floating island requires bridging/support, which is unavailable in nonplanar mode");
            for (int wall = 0; wall < settings.walls; ++wall) {
                const auto contours = inset(safe,wall*settings.width,tol);
                for (const auto &p : contours) if (p.size() >= 3)
                    layer.paths.push_back(planner.path(p,true,k,wall ? Role::InnerWall : Role::OuterWall));
            }
            Paths64 interior = inset(safe,settings.walls*settings.width,tol);
            // Skin propagation follows material presence on neighboring supplied
            // bands, independently for each object and all its islands/holes.
            Paths64 top = interior, bottom = interior;
            for (int n = 1; n <= settings.top_layers; ++n) {
                if (k+size_t(n) >= bands.size()) { top.clear(); break; }
                top = intersect(top,bands[k+n]);
            }
            top = settings.top_layers ? difference(interior,top) : Paths64{};
            for (int n = 1; n <= settings.bottom_layers; ++n) {
                if (size_t(n) > k) { bottom.clear(); break; }
                bottom = intersect(bottom,bands[k-n]);
            }
            bottom = settings.bottom_layers ? difference(interior,bottom) : Paths64{};
            const Paths64 bottom_only = difference(bottom,top), solid = Union(top,bottom,FillRule::NonZero);
            const Paths64 sparse = difference(interior,solid);
            auto fill = [&](const Paths64 &region, double spacing, Role role) {
                for (const auto &p : infill(region,stack,spacing,settings.infill_angle+90.*(k%2))) if (p.size() >= 2)
                    layer.paths.push_back(planner.path(p,false,k,role));
            };
            fill(top,settings.width,Role::TopSkin);
            fill(bottom_only,settings.width,Role::BottomSkin);
            if (settings.density > 0) fill(sparse,settings.width/settings.density,
                                          settings.density == 1. ? Role::SolidInfill : Role::SparseInfill);
        }
    }
    if (std::none_of(output.begin(),output.end(),[](const Layer &l) { return !l.paths.empty(); }))
        invalid("no printable bead fits in the supplied material bands");
    for (Layer &layer : output) if (!layer.paths.empty()) {
        layer.min_z = std::numeric_limits<double>::infinity(); layer.max_z = -layer.min_z;
        auto add = [&](XYZ p) { layer.min_z = std::min(layer.min_z,p.z); layer.max_z = std::max(layer.max_z,p.z); };
        for (const auto &path : layer.paths) {
            add(path.start);
            for (const auto &segment : path.segments) add(segment.end);
        }
        if (layer.min_z < 0. || layer.max_z > settings.max_z) invalid("generated path exceeds machine Z limits");
    }
    return output;
}

} // namespace Slic3r::Nonplanar
