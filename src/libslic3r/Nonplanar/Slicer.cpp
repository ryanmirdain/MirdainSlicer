#include "Slicer.hpp"
#include "../I18N.hpp"
#include "../format.hpp"
#include <clipper2/clipper.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
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
                double lo = t.vertices.front().z, hi = lo;
                for (XYZ p : t.vertices) { lo = std::min(lo,p.z); hi = std::max(hi,p.z); }
                // The surface facets already lie within this cell. Intersect
                // the original triangle, avoiding extra clipping seams.
                cells[j*(stack.x.size()-1)+i].push_back({std::vector<XYZ>(t.vertices.begin(),t.vertices.end()),normal,lo,hi});
            }
        }
        for (auto &cell : cells)
            std::sort(cell.begin(),cell.end(),[](const Piece &a, const Piece &b) { return a.max_z < b.max_z; });
    }

    Paths64 cross_section(const Surface &surface)
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
            patches[cell] = stack.patch(surface,i,j);
            for (double &height : patches[cell].control)
                height -= boundary_epsilon;
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
        return loops;
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
        if (enter < leave) {
            XYZ a = mix(hits[0], hits[1], enter), b = mix(hits[0], hits[1], leave);
            const XYZ surface_normal = cross(subtract(triangle[1], triangle[0]), subtract(triangle[2], triangle[0]));
            const XYZ direction = cross(surface_normal, piece.normal);
            if (dot(subtract(b, a), direction) < 0.) std::swap(a, b);
            emit(xy(a), xy(b));
        }
    }
};

} // namespace

Surface Job::surface_at(double z) const
{
    const size_t last = stack.layers.size() - 1;
    const Surface *below, *above;
    double fraction, extension = 0.;
    if (z <= 0.) {
        below = above = &stack.bed;
        fraction = 0.; extension = z;
    } else if (z >= layer_z(last)) {
        below = above = &stack.layers[last];
        fraction = 0.; extension = z - layer_z(last);
    } else {
        const size_t layer = size_t(std::max(0., std::ceil((z - first_layer_height) / layer_height)));
        below = layer ? &stack.layers[layer - 1] : &stack.bed;
        above = &stack.layers[layer];
        fraction = (z - (layer ? layer_z(layer - 1) : 0.)) / (layer ? layer_height : first_layer_height);
    }
    Surface surface;
    surface.plane = false;
    surface.patches.reserve(below->patches.size());
    for (size_t cell = 0; cell < below->patches.size(); ++cell) {
        Patch patch = below->patches[cell];
        for (size_t n = 0; n < patch.control.size(); ++n)
            patch.control[n] += fraction * (above->patches[cell].control[n] - patch.control[n]) + extension;
        surface.patches.push_back(std::move(patch));
    }
    return surface;
}

Vec3d Job::to_model(const Vec3d &point) const
{
    // Match the historical height-profile inverse: clamp the XY lookup only,
    // retain the requested XY position and extend vertically at unit rate.
    const XY xy {std::clamp(point.x(), stack.x.front(), stack.x.back()),
                 std::clamp(point.y(), stack.y.front(), stack.y.back())};
    const size_t last = stack.layers.size() - 1;
    double z;
    if (point.z() <= 0.) {
        z = stack.value(stack.bed, xy) + point.z();
    } else if (point.z() >= layer_z(last)) {
        z = stack.value(stack.layers[last], xy) + point.z() - layer_z(last);
    } else {
        const size_t layer = size_t(std::max(0., std::ceil((point.z() - first_layer_height) / layer_height)));
        const double fraction = (point.z() - (layer ? layer_z(layer - 1) : 0.)) /
                                (layer ? layer_height : first_layer_height);
        const double below = stack.value(layer ? stack.layers[layer - 1] : stack.bed, xy);
        z = below + fraction * (stack.value(stack.layers[layer], xy) - below);
    }
    return {point.x(), point.y(), z};
}

std::vector<Polygons> cross_sections(const Job &job, const Object &object, const std::vector<float> &zs,
                                    const std::function<void()> &cancel)
{
    const auto &stack = job.stack;
    for (const auto &triangle : object.triangles)
        for (const auto &p : triangle.vertices)
            if (!stack.covers(xy(p)))
                invalid(_u8L("shared grid does not cover the placed model footprint"));
    Intersection intersection(stack, object, cancel);
    std::vector<Polygons> result;
    result.reserve(zs.size());
    for (float z : zs) {
        cancel();
        Polygons polygons;
        for (const auto &contour : intersection.cross_section(job.surface_at(z))) {
            Points points;
            for (const auto &p : contour)
                points.emplace_back(p.x, p.y);
            polygons.emplace_back(std::move(points));
        }
        result.push_back(std::move(polygons));
    }
    return result;
}
} // namespace Slic3r::Nonplanar
