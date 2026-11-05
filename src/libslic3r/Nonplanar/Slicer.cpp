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

// A mesh face intersected with a bilinear height patch has the XY equation
// A + B*u + C*v + D*u*v = 0. Trace that curve, not triangle subdivisions.
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
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            cancel();
            if (cells[cell].empty()) continue;
            const size_t i = cell%nx, j = cell/nx;
            const Patch lower = stack.patch(below,i,j), upper = stack.patch(above,i,j);
            const Patch h {lower.x,lower.y,(lower.c+upper.c)*.5-boundary_epsilon,
                           (lower.ax+upper.ax)*.5,(lower.ay+upper.ay)*.5,(lower.axy+upper.axy)*.5};
            const std::array<double,4> heights {h.value({stack.x[i],stack.y[j]}),h.value({stack.x[i+1],stack.y[j]}),
                h.value({stack.x[i],stack.y[j+1]}),h.value({stack.x[i+1],stack.y[j+1]})};
            const double low = *std::min_element(heights.begin(),heights.end());
            const double high = *std::max_element(heights.begin(),heights.end());
            auto first = std::lower_bound(cells[cell].begin(),cells[cell].end(),low,
                [](const Piece &p, double z) { return p.max_z < z; });
            for (auto it = first; it != cells[cell].end(); ++it) {
                const Piece &piece = *it;
                if (piece.min_z > high) continue;
                trace(piece,h,line);
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
    template<class Emit> void trace(const Piece &piece, const Patch &h, Emit emit)
    {
        std::vector<XYZ> hits;
        auto add = [&](XYZ p) {
            for (XYZ q : hits) if (distance(p,q) < 1e-8) return;
            hits.push_back(p);
        };
        for (size_t n = 0; n < piece.polygon.size(); ++n) {
            const XYZ a = piece.polygon[n], b = piece.polygon[(n+1)%piece.polygon.size()];
            auto f = [&](double t) { XYZ p = mix(a,b,t); return p.z-h.value(xy(p)); };
            const double c = f(0.), end = f(1.), mid = f(.5);
            const double qa = 2.*(end+c-2.*mid), qb = end-c-qa;
            auto root = [&](double t) {
                if (t >= -1e-10 && t <= 1.+1e-10) add(mix(a,b,std::clamp(t,0.,1.)));
            };
            if (std::abs(c) < 1e-10) add(a);
            if (std::abs(end) < 1e-10) add(b);
            if (std::abs(qa) < 1e-13) {
                if (std::abs(qb) > 1e-13) root(-c/qb);
            } else {
                const double discriminant = qb*qb-4.*qa*c;
                if (discriminant >= 0.) {
                    const double q = -.5*(qb+std::copysign(std::sqrt(discriminant),qb));
                    if (q != 0.) { root(q/qa); root(c/q); }
                    else root(-qb/(2.*qa));
                }
            }
        }
        if (hits.size() < 2) return;
        const XYZ n = piece.normal, origin = piece.polygon.front();
        const double A = n.x*(h.x-origin.x)+n.y*(h.y-origin.y)+n.z*(h.c-origin.z);
        const double B = n.x+n.z*h.ax, C = n.y+n.z*h.ay, D = n.z*h.axy;
        // Choose the graph direction with the greater boundary extent. Both
        // directions describe the same conic; vertical faces remain supported.
        double xmin = hits[0].x, xmax = xmin, ymin = hits[0].y, ymax = ymin;
        for (XYZ p : hits) {
            xmin = std::min(xmin,p.x); xmax = std::max(xmax,p.x);
            ymin = std::min(ymin,p.y); ymax = std::max(ymax,p.y);
        }
        const bool along_x = xmax-xmin >= ymax-ymin;
        const double linear = along_x ? B : C, denominator = along_x ? C : B;
        auto coordinate = [&](XYZ p) { return along_x ? p.x-h.x : p.y-h.y; };
        auto at = [&](double t) {
            const double other = -(A+linear*t)/(denominator+D*t);
            return along_x ? XY{h.x+t,h.y+other} : XY{h.x+other,h.y+t};
        };
        auto inside = [&](XY p) {
            const XYZ q {p.x,p.y,h.value(p)};
            for (size_t v = 0; v < piece.polygon.size(); ++v) {
                const XYZ a = piece.polygon[v], b = piece.polygon[(v+1)%piece.polygon.size()];
                if (dot(cross(subtract(b,a),subtract(q,a)),n) < -1e-8) return false;
            }
            return true;
        };
        // A factored conic consists of two straight branches. Handle each
        // branch separately so the graph's removable pole does not lose a cut.
        if (D != 0. && std::abs(A*D-B*C) <= 1e-14*(std::abs(A*D)+std::abs(B*C)+1e-20)) {
            for (int axis = 0; axis < 2; ++axis) {
                const double fixed = axis ? -B/D+h.y : -C/D+h.x;
                std::vector<XYZ> branch;
                for (XYZ p : hits)
                    if (std::abs((axis ? p.y : p.x)-fixed) < 1e-7) branch.push_back(p);
                std::sort(branch.begin(),branch.end(),[&](XYZ a, XYZ b) { return axis ? a.x < b.x : a.y < b.y; });
                for (size_t v = 1; v < branch.size(); ++v) {
                    const XY a = xy(branch[v-1]), b = xy(branch[v]);
                    if (inside({(a.x+b.x)*.5,(a.y+b.y)*.5})) emit(a,b);
                }
            }
            return;
        }
        std::sort(hits.begin(),hits.end(),[&](XYZ a, XYZ b) { return coordinate(a) < coordinate(b); });
        for (size_t v = 1; v < hits.size(); ++v) {
            const double a = coordinate(hits[v-1]), b = coordinate(hits[v]);
            if (b-a < 1e-10) continue;
            if (D != 0. && -denominator/D > a && -denominator/D < b) continue;
            const XY mid = at((a+b)*.5);
            if (!std::isfinite(mid.x) || !std::isfinite(mid.y) || !inside(mid)) continue;
            // Sample only the intersection curve. For this rational graph the
            // maximum chord deviation occurs where its derivative matches the
            // chord slope; evaluate that point rather than subdividing faces.
            auto curve = [&](auto &&self, double ta, double tb, XY pa, XY pb, unsigned depth) -> void {
                double split = (ta+tb)*.5;
                const double da = denominator+D*ta, db = denominator+D*tb;
                if (D != 0. && da*db > 0.)
                    split = (std::copysign(std::sqrt(da*db),da)-denominator)/D;
                const XY pm = at(split);
                const double dx = pb.x-pa.x, dy = pb.y-pa.y;
                const double length = std::hypot(dx,dy);
                const double error = length ? std::abs(dx*(pm.y-pa.y)-dy*(pm.x-pa.x))/length : 0.;
                if (error > stack.limits.geometry_tolerance*.25 && depth < 32) {
                    self(self,ta,split,pa,pm,depth+1); self(self,split,tb,pm,pb,depth+1);
                } else emit(pa,pb);
            };
            curve(curve,a,b,xy(hits[v-1]),xy(hits[v]),0);
        }
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

    double account_volume(XY a, XY b, size_t k, double &projected_area) const
    {
        const double length = std::hypot(b.x-a.x,b.y-a.y);
        const XY mid {(a.x+b.x)*.5,(a.y+b.y)*.5};
        // Within a grid cell the vertical gap along a straight XY segment is
        // quadratic. Simpson integration is exact, with no bead-footprint booleans.
        const double gap = (stack.vertical_gap(k,a)+4.*stack.vertical_gap(k,mid)+stack.vertical_gap(k,b))/6.;
        projected_area = settings.width*length;
        return projected_area*gap;
    }

    void segment(Path &out, XY a, XY b, size_t k, unsigned depth)
    {
        const XY mid {(a.x+b.x)*.5,(a.y+b.y)*.5};
        const XYZ p {a.x,a.y,stack.value(stack.layers[k],a)}, q {b.x,b.y,stack.value(stack.layers[k],b)};
        const double zm = stack.value(stack.layers[k],mid), gap = stack.vertical_gap(k,mid);
        const double chord_error = std::abs(zm-(p.z+q.z)*.5);
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
