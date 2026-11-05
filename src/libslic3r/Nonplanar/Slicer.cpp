#include "Slicer.hpp"
#include "../I18N.hpp"
#include "../format.hpp"
#include <clipper2/clipper.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace Slic3r::Nonplanar {
namespace {
using namespace Clipper2Lib;
constexpr double scale = 1e6;
constexpr double boundary_epsilon = 1e-8;
[[noreturn]] void invalid(const std::string &s) { throw std::invalid_argument(Slic3r::format(_u8L("Nonplanar slicing: %1%"),s)); }
XY xy(XYZ p) { return {p.x,p.y}; }
XYZ mix(XYZ a, XYZ b, double t)
{ return {a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t,a.z+(b.z-a.z)*t}; }
double distance(XYZ a, XYZ b)
{ return std::hypot(std::hypot(a.x-b.x,a.y-b.y),a.z-b.z); }
Point64 point(XY p) { return {std::llround(p.x*scale),std::llround(p.y*scale)}; }
XY real(Point64 p) { return {p.x/scale,p.y/scale}; }
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

// Intersect the original mesh with a tolerance-bounded tessellation of the
// shared spline. Only the surface is approximated; model triangles stay intact.
class Intersection {
    struct Piece {
        std::vector<XYZ> polygon;
        XYZ normal;
        double min_z, max_z;
    };
    const LayerStack &stack;
    const std::function<void()> &cancel;
    std::vector<std::vector<Piece>> cells;

    static XYZ cross(XYZ a, XYZ b)
    { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
    static XYZ subtract(XYZ a, XYZ b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
    static double dot(XYZ a, XYZ b) { return a.x*b.x+a.y*b.y+a.z*b.z; }

public:
    Intersection(const LayerStack &s, const Object &object, const std::function<void()> &c)
        : stack(s), cancel(c), cells((s.x.size()-1)*(s.y.size()-1))
    {
        for (const Triangle &t : object.triangles) {
            cancel();
            XYZ normal = cross(subtract(t.vertices[1],t.vertices[0]),subtract(t.vertices[2],t.vertices[0]));
            const double length = std::sqrt(dot(normal,normal));
            if (length < 1e-16) continue;
            normal = {normal.x/length,normal.y/length,normal.z/length};
            double xmin = t.vertices[0].x, xmax = xmin, ymin = t.vertices[0].y, ymax = ymin;
            for (XYZ p : t.vertices) {
                xmin = std::min(xmin,p.x); xmax = std::max(xmax,p.x);
                ymin = std::min(ymin,p.y); ymax = std::max(ymax,p.y);
            }
            const auto [i0,j0] = stack.cell({xmin,ymin});
            const auto [i1,j1] = stack.cell({xmax,ymax});
            for (size_t j = j0; j <= j1; ++j) for (size_t i = i0; i <= i1; ++i) {
                std::vector<XYZ> polygon(t.vertices.begin(),t.vertices.end());
                polygon = clip_axis(polygon,0,stack.x[i],true); polygon = clip_axis(polygon,0,stack.x[i+1],false);
                polygon = clip_axis(polygon,1,stack.y[j],true); polygon = clip_axis(polygon,1,stack.y[j+1],false);
                if (polygon.size() < 3) continue;
                XYZ sum {0.,0.,0.};
                for (size_t v = 1; v+1 < polygon.size(); ++v) {
                    const XYZ c = cross(subtract(polygon[v],polygon[0]),subtract(polygon[v+1],polygon[0]));
                    sum = {sum.x+c.x,sum.y+c.y,sum.z+c.z};
                }
                if (dot(sum,sum) < 1e-24) continue;
                double lo = polygon.front().z, hi = lo;
                for (XYZ p : polygon) { lo = std::min(lo,p.z); hi = std::max(hi,p.z); }
                cells[j*(stack.x.size()-1)+i].push_back({std::move(polygon),normal,lo,hi});
            }
        }
        for (auto &cell : cells)
            std::sort(cell.begin(),cell.end(),[](const Piece &a, const Piece &b) { return a.max_z < b.max_z; });
    }

    Paths64 cross_section(const Surface &below, const Surface &above)
    {
        // Join curve segments just as a planar mesh slicer joins triangle cuts.
        struct Node { XY p; std::vector<size_t> edges; };
        struct Edge { size_t a,b; bool used = false; };
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::map<std::pair<int64_t,int64_t>,std::vector<size_t>> buckets;
        std::set<std::pair<size_t,size_t>> unique;
        constexpr double join = 2e-6;
        auto node = [&](XY p) {
            const int64_t x = std::llround(p.x/join), y = std::llround(p.y/join);
            for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) {
                auto it = buckets.find({x+dx,y+dy});
                if (it != buckets.end()) for (size_t n : it->second)
                    if (std::hypot(nodes[n].p.x-p.x,nodes[n].p.y-p.y) <= join) return n;
            }
            const size_t n = nodes.size();
            nodes.push_back({p,{}}); buckets[{x,y}].push_back(n); return n;
        };
        auto line = [&](XY a, XY b) {
            const size_t na = node(a), nb = node(b);
            if (na == nb || !unique.insert(std::minmax(na,nb)).second) return;
            nodes[na].edges.push_back(edges.size()); nodes[nb].edges.push_back(edges.size());
            edges.push_back({na,nb});
        };
        const size_t nx = stack.x.size()-1;
        std::vector<Patch> patches(cells.size());
        std::vector<size_t> steps_x(nx,1), steps_y(stack.y.size()-1,1);
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            if (cells[cell].empty()) continue;
            const size_t i = cell%nx, j = cell/nx;
            patches[cell] = stack.patch(below,i,j);
            const Patch upper = stack.patch(above,i,j);
            for (size_t n = 0; n < 16; ++n)
                patches[cell].control[n] = .5*(patches[cell].control[n]+upper.control[n])-boundary_epsilon;
            const size_t steps = std::max(size_t(1),size_t(std::ceil(std::sqrt(patches[cell].linear_error()/(stack.limits.geometry_tolerance*.25)))));
            steps_x[i] = std::max(steps_x[i],steps);
            steps_y[j] = std::max(steps_y[j],steps);
        }
        // Neighbors share edge subdivision counts along each grid axis, so
        // contours have no hanging vertices. Hessian bounds control height
        // error without recursively splitting the original mesh faces.
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            cancel();
            if (cells[cell].empty()) continue;
            const Patch &h = patches[cell];
            const auto [low,high] = h.bounds();
            auto first = std::lower_bound(cells[cell].begin(),cells[cell].end(),low,
                [](const Piece &p, double z) { return p.max_z < z; });
            if (first == cells[cell].end() || std::none_of(first,cells[cell].end(),[&](const Piece &p) { return p.min_z <= high; })) continue;
            const size_t sx = steps_x[cell%nx], sy = steps_y[cell/nx];
            std::vector<XYZ> vertices((sx+1)*(sy+1));
            for (size_t v = 0; v <= sy; ++v) for (size_t u = 0; u <= sx; ++u) {
                const XY p {h.x+h.dx*double(u)/sx,h.y+h.dy*double(v)/sy};
                vertices[v*(sx+1)+u] = {p.x,p.y,h.value(p)};
            }
            auto facet = [&](const std::array<XYZ,3> &triangle) {
                const double lo = std::min({triangle[0].z,triangle[1].z,triangle[2].z});
                const double hi = std::max({triangle[0].z,triangle[1].z,triangle[2].z});
                const auto begin = std::lower_bound(first,cells[cell].end(),lo,
                    [](const Piece &p, double z) { return p.max_z < z; });
                for (auto it = begin; it != cells[cell].end(); ++it)
                    if (it->min_z <= hi) trace(*it,triangle,line);
            };
            for (size_t v = 0; v < sy; ++v) for (size_t u = 0; u < sx; ++u) {
                const size_t n = v*(sx+1)+u;
                facet({vertices[n],vertices[n+1],vertices[n+sx+2]});
                facet({vertices[n],vertices[n+sx+2],vertices[n+sx+1]});
            }
        }
        Paths64 loops;
        for (size_t e = 0; e < edges.size(); ++e) {
            if (edges[e].used) continue;
            const size_t start = edges[e].a;
            size_t current = start, next_edge = e;
            Path64 path;
            do {
                path.push_back(point(nodes[current].p));
                Edge &edge = edges[next_edge]; edge.used = true;
                current = edge.a == current ? edge.b : edge.a;
                if (current == start) break;
                auto next = std::find_if(nodes[current].edges.begin(),nodes[current].edges.end(),
                    [&](size_t n) { return !edges[n].used; });
                if (next == nodes[current].edges.end())
                    invalid(_u8L("mesh intersection contains an open contour"));
                next_edge = *next;
            } while (current != start);
            if (path.size() >= 3) loops.push_back(std::move(path));
        }
        return Union(loops,FillRule::EvenOdd);
    }

private:
    template<class Emit> void trace(const Piece &piece, const std::array<XYZ,3> &triangle, Emit emit)
    {
        const XYZ origin = piece.polygon.front();
        std::array<double,3> residual;
        for (size_t n = 0; n < 3; ++n) residual[n] = dot(subtract(triangle[n],origin),piece.normal);
        if (*std::min_element(residual.begin(),residual.end()) > 0. ||
            *std::max_element(residual.begin(),residual.end()) < 0.) return;
        std::vector<XYZ> hits;
        auto add = [&](XYZ p) {
            for (XYZ q : hits) if (distance(p,q) < 1e-10) return;
            hits.push_back(p);
        };
        for (size_t n = 0; n < 3; ++n) {
            const size_t next = (n+1)%3;
            if (std::abs(residual[n]) < 1e-12) add(triangle[n]);
            if ((residual[n] >= 0.) != (residual[next] >= 0.))
                add(mix(triangle[n],triangle[next],residual[n]/(residual[n]-residual[next])));
        }
        if (hits.size() != 2) return;
        double enter = 0., leave = 1.;
        for (size_t n = 0; n < piece.polygon.size(); ++n) {
            const XYZ a = piece.polygon[n], edge = subtract(piece.polygon[(n+1)%piece.polygon.size()],a);
            const double f0 = dot(cross(edge,subtract(hits[0],a)),piece.normal);
            const double f1 = dot(cross(edge,subtract(hits[1],a)),piece.normal);
            if (f0 < 0. && f1 < 0.) return;
            if ((f0 >= 0.) != (f1 >= 0.)) {
                const double t = f0/(f0-f1);
                if (f0 < 0.) enter = std::max(enter,t); else leave = std::min(leave,t);
            }
        }
        if (enter < leave) emit(xy(mix(hits[0],hits[1],enter)),xy(mix(hits[0],hits[1],leave)));
    }
};

class Planner {
public:
    Planner(const LayerStack &s, const Settings &p) : stack(s), settings(p) {}

    Path path(const Path64 &points, bool closed, size_t k, Role role, uint64_t object_id)
    {
        Path out;
        out.role = role;
        out.object_id = object_id;
        const size_t seam = closed ? size_t(std::min_element(points.begin(),points.end(),[](Point64 a, Point64 b) {
            return a.x == b.x ? a.y < b.y : a.x < b.x;
        })-points.begin()) : 0;
        const XY first = real(points[seam]);
        out.start = {first.x,first.y,stack.value(stack.layers[k],first)};
        for (size_t n = 1; n < points.size()+size_t(closed); ++n) {
            const XY a = real(points[(seam+n-1)%points.size()]), b = real(points[(seam+n)%points.size()]);
            // Split at every grid seam before curvature/flow subdivision. A
            // spline span is polynomial within a cell, never across spans.
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

    double account_volume(XY a, XY b, size_t k, double &projected_area) const
    {
        const double length = std::hypot(b.x-a.x,b.y-a.y);
        // A bicubic span restricted to a line has degree six. Four-point
        // Gauss-Legendre integrates its gap exactly between grid seams.
        constexpr std::array<double,4> nodes {.06943184420297371,.33000947820757187,.6699905217924281,.9305681557970262};
        constexpr std::array<double,4> weights {.17392742256872693,.32607257743127307,.32607257743127307,.17392742256872693};
        double gap = 0.;
        for (size_t n = 0; n < 4; ++n)
            gap += weights[n]*stack.vertical_gap(k,{a.x+(b.x-a.x)*nodes[n],a.y+(b.y-a.y)*nodes[n]});
        projected_area = settings.width*length;
        return projected_area*gap;
    }

    void segment(Path &out, XY a, XY b, size_t k, unsigned depth)
    {
        const XY mid {(a.x+b.x)*.5,(a.y+b.y)*.5};
        const XYZ p {a.x,a.y,stack.value(stack.layers[k],a)}, q {b.x,b.y,stack.value(stack.layers[k],b)};
        const double zm = stack.value(stack.layers[k],mid), gap = stack.vertical_gap(k,mid);
        const auto [i,j] = stack.cell(mid);
        const auto curve = stack.patch(stack.layers[k],i,j).line(a,b);
        double chord_error = 0.;
        for (size_t n = 0; n < curve.size(); ++n)
            chord_error = std::max(chord_error,std::abs(curve[n]-(p.z+(q.z-p.z)*double(n)/6.)));
        if (depth < 32 && (distance(p,q) > stack.limits.max_segment_length-.00018 ||
                           chord_error > stack.limits.chord_tolerance*.5)) {
            segment(out,a,mid,k,depth+1); segment(out,mid,b,k,depth+1); return;
        }
        const double ds = std::hypot(b.x-a.x,b.y-a.y);
        if (ds < 1e-10) return;
        const double normal = stack.normal_spacing(k,mid);
        const XY grad = stack.gradient(stack.layers[k],mid);
        const double cross_slope = (-grad.x*(b.y-a.y)+grad.y*(b.x-a.x))/ds;
        // Preview uses the local tangent band; cross-width approximation is
        // independent of whether the centerline can be sliced and emitted.
        const double max_top_z = std::max({p.z,q.z,zm})+settings.width*.5*std::abs(cross_slope);
        double projected_area;
        const double volume = account_volume(a,b,k,projected_area);
        out.projected_area += projected_area;
        out.segments.push_back({q,volume,gap,normal,settings.width*std::sqrt(1.+cross_slope*cross_slope),-cross_slope,max_top_z});
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
    if (settings.walls < 0 || settings.width <= 0 || settings.density < 0 || settings.density > 1)
        invalid(_u8L("a nonnegative wall count, a positive line width and density in [0,1] are required"));
    const double tol = stack.limits.geometry_tolerance;
    for (const Object &object : objects) for (const Triangle &t : object.triangles) for (XYZ p : t.vertices) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) invalid(_u8L("model contains nonfinite coordinates"));
        if (!stack.covers(xy(p))) invalid(_u8L("shared grid does not cover the placed model footprint"));
    }
    std::vector<Layer> output(stack.layers.size());
    Planner planner(stack,settings);
    for (const Object &object : objects) {
        Intersection engine(stack,object,cancel);
        std::vector<Paths64> sections;
        for (size_t k = 0; k < stack.layers.size(); ++k) {
            cancel();
            const Surface &below = k ? stack.layers[k-1] : stack.bed;
            sections.push_back(engine.cross_section(below,stack.layers[k]));
        }
        for (size_t k = 0; k < sections.size(); ++k) {
            cancel();
            Layer &layer = output[k];
            const double margin = settings.width*.5+tol;
            Paths64 safe = inset(sections[k],margin,tol);
            for (int wall = 0; wall < settings.walls; ++wall) {
                const auto contours = inset(safe,wall*settings.width,tol);
                for (const auto &p : contours) if (p.size() >= 3)
                    layer.paths.push_back(planner.path(p,true,k,wall ? Role::InnerWall : Role::OuterWall,object.id));
            }
            Paths64 interior = inset(safe,std::max(0.,settings.walls-.5)*settings.width,tol);
            // Skin propagation follows material presence on neighboring supplied
            // bands, independently for each object and all its islands/holes.
            Paths64 top = interior, bottom = interior;
            for (int n = 1; n <= settings.top_layers; ++n) {
                if (k+size_t(n) >= sections.size()) { top.clear(); break; }
                top = intersect(top,sections[k+n]);
            }
            top = settings.top_layers ? difference(interior,top) : Paths64{};
            for (int n = 1; n <= settings.bottom_layers; ++n) {
                if (size_t(n) > k) { bottom.clear(); break; }
                bottom = intersect(bottom,sections[k-n]);
            }
            bottom = settings.bottom_layers ? difference(interior,bottom) : Paths64{};
            const Paths64 bottom_only = difference(bottom,top), solid = Union(top,bottom,FillRule::NonZero);
            const Paths64 sparse = difference(interior,solid);
            auto fill = [&](const Paths64 &region, double spacing, Role role) {
                // Flat-ended strip beads must overlap the innermost wall and
                // the skin partition boundary. Dilation includes every grid
                // line whose width covers the target; clipping preserves the
                // sampled cross-section boundary and its holes.
                const auto expanded = InflatePaths(region,settings.width*.5*scale,Clipper2Lib::JoinType::Round,
                                                   Clipper2Lib::EndType::Polygon,2.,tol*scale/8.);
                const auto centers = intersect(expanded,safe);
                for (const auto &p : infill(centers,stack,spacing,settings.infill_angle+90.*(k%2))) if (p.size() >= 2)
                    layer.paths.push_back(planner.path(p,false,k,role,object.id));
            };
            fill(top,settings.width,Role::TopSkin);
            fill(bottom_only,settings.width,Role::BottomSkin);
            if (settings.density > 0) fill(sparse,settings.width/settings.density,
                                          settings.density == 1. ? Role::SolidInfill : Role::SparseInfill);
        }
    }
    if (std::none_of(output.begin(),output.end(),[](const Layer &l) { return !l.paths.empty(); }))
        invalid(_u8L("no printable bead fits in the supplied material bands"));
    for (Layer &layer : output) if (!layer.paths.empty()) {
        layer.min_z = std::numeric_limits<double>::infinity(); layer.max_z = -layer.min_z;
        auto add = [&](XYZ p) { layer.min_z = std::min(layer.min_z,p.z); layer.max_z = std::max(layer.max_z,p.z); };
        for (const auto &path : layer.paths) {
            add(path.start);
            for (const auto &segment : path.segments) {
                add(segment.end);
                // Machine bounds apply to emitted coordinates, not preview envelopes.
            }
        }
        if (layer.min_z < 0. || layer.max_z > settings.max_z) invalid(_u8L("generated path exceeds machine Z limits"));
    }
    return output;
}

} // namespace Slic3r::Nonplanar
