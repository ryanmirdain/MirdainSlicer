#include "Slicer.hpp"
#include "../I18N.hpp"
#include "../format.hpp"
#include <clipper2/clipper.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <chrono>
#include <memory>
#include <mutex>
#include <optional>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

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
        std::array<XYZ,3> polygon;
        std::array<size_t,3> vertex_ids;
        size_t face_id;
        XYZ normal;
        double min_z, max_z;
    };
    const LayerStack &stack;
    std::vector<std::vector<Piece>> cells;
    // Pieces are sorted by max_z, so a facet's candidates are a suffix whose
    // min_z must still be tested. Contiguous min_z values with per-block
    // minima skip pieces that start above the facet without changing the
    // visiting order, which fixes contour node order and output.
    static constexpr size_t z_block = 16;
    struct ZIndex { std::vector<double> min_z, block_min_z; };
    std::vector<ZIndex> z_index;
    using MeshEdge = std::pair<size_t,size_t>;
    std::map<MeshEdge,MeshEdge> collapsed_edges;

    MeshEdge edge_root(MeshEdge edge) const
    {
        auto it = collapsed_edges.find(edge);
        while (it != collapsed_edges.end()) {
            edge = it->second;
            it = collapsed_edges.find(edge);
        }
        return edge;
    }

    static XYZ cross(XYZ a, XYZ b)
    { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
    static XYZ subtract(XYZ a, XYZ b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
    static double dot(XYZ a, XYZ b) { return a.x*b.x+a.y*b.y+a.z*b.z; }

public:
    Intersection(const LayerStack &s, const Object &object, const std::function<void()> &c)
        : stack(s), cells((s.x.size()-1)*(s.y.size()-1))
    {
        for (size_t face = 0; face < object.triangles.size(); ++face) {
            const Triangle &t = object.triangles[face];
            c();
            XYZ normal = cross(subtract(t.vertices[1],t.vertices[0]),subtract(t.vertices[2],t.vertices[0]));
            const double length = std::sqrt(dot(normal,normal));
            if (length < 1e-16) {
                // A collinear face has no area to intersect, but its edges
                // still connect the neighboring faces along the same line.
                std::array<MeshEdge,3> roots;
                for (size_t n = 0; n < 3; ++n) {
                    const auto ids = std::minmax(t.vertex_ids[n],t.vertex_ids[(n+1)%3]);
                    roots[n] = edge_root({ids.first,ids.second});
                }
                const MeshEdge root = *std::min_element(roots.begin(),roots.end());
                for (MeshEdge edge : roots) if (edge != root) collapsed_edges[edge] = root;
                continue;
            }
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
                cells[j*(stack.x.size()-1)+i].push_back({t.vertices,t.vertex_ids,face,normal,lo,hi});
            }
        }
        z_index.resize(cells.size());
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            auto &pieces = cells[cell];
            std::sort(pieces.begin(),pieces.end(),[](const Piece &a, const Piece &b) { return a.max_z < b.max_z; });
            auto &index = z_index[cell];
            index.min_z.reserve(pieces.size());
            for (const Piece &piece : pieces) index.min_z.push_back(piece.min_z);
            for (size_t n = 0; n < pieces.size(); n += z_block)
                index.block_min_z.push_back(*std::min_element(index.min_z.begin()+n,
                    index.min_z.begin()+std::min(n+z_block,pieces.size())));
        }
    }

    struct Tessellation {
        struct Cell {
            double low = 0., high = 0.;
            size_t sx = 0, sy = 0;
            std::vector<XYZ> vertices;
        };
        std::vector<Cell> cells;
    };

    Tessellation tessellate(const Surface &surface, double min_offset, double max_offset,
                            const std::function<void()> &cancel) const
    {
        const size_t nx = stack.x.size()-1;
        std::vector<Patch> patches(cells.size());
        std::vector<size_t> steps_x(nx,1), steps_y(stack.y.size()-1,1);
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            if (cells[cell].empty()) continue;
            const size_t i = cell%nx, j = cell/nx;
            patches[cell] = stack.patch(surface,i,j);
            for (double &height : patches[cell].control) height -= boundary_epsilon;
            const size_t steps = std::max(size_t(1),size_t(std::ceil(std::sqrt(patches[cell].linear_error()/(stack.limits.geometry_tolerance*.25)))));
            steps_x[i] = std::max(steps_x[i],steps);
            steps_y[j] = std::max(steps_y[j],steps);
        }
        Tessellation result;
        result.cells.resize(cells.size());
        // Conforming XY subdivision is invariant under vertical translation.
        // Evaluate it once, then reuse the vertices for all three contours.
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            cancel();
            if (cells[cell].empty()) continue;
            const Patch &h = patches[cell];
            auto &out = result.cells[cell];
            std::tie(out.low,out.high) = h.bounds();
            auto first = std::lower_bound(cells[cell].begin(),cells[cell].end(),out.low+min_offset,
                [](const Piece &p, double z) { return p.max_z < z; });
            if (first == cells[cell].end() || std::none_of(first,cells[cell].end(),
                [high = out.high+max_offset](const Piece &p) { return p.min_z <= high; })) continue;
            out.sx = steps_x[cell%nx]; out.sy = steps_y[cell/nx];
            out.vertices.resize((out.sx+1)*(out.sy+1));
            for (size_t v = 0; v <= out.sy; ++v) for (size_t u = 0; u <= out.sx; ++u) {
                const XY p {h.x+h.dx*double(u)/out.sx,h.y+h.dy*double(v)/out.sy};
                out.vertices[v*(out.sx+1)+u] = {p.x,p.y,h.value(p)};
            }
        }
        return result;
    }

    Paths64 cross_section(const Tessellation &surface, double offset, const std::function<void()> &cancel) const
    {
        // Coordinates reconcile roundoff only within the same indexed mesh
        // feature. Coincident boundaries of different shells remain distinct.
        struct Node { XY p; std::vector<size_t> edges; };
        struct Edge { size_t a,b; bool used = false; };
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        struct BucketKey {
            std::array<size_t,3> feature; int64_t x, y;
            bool operator==(const BucketKey &o) const { return feature == o.feature && x == o.x && y == o.y; }
        };
        struct BucketHash {
            size_t operator()(const BucketKey &k) const {
                size_t h = std::hash<int64_t>()(k.x)*0x9E3779B97F4A7C15ull ^ std::hash<int64_t>()(k.y);
                for (size_t v : k.feature) h = (h ^ v)*0x100000001B3ull;
                return h;
            }
        };
        // Lookups only; never iterated, so hashing cannot change node order.
        std::unordered_map<BucketKey,std::vector<size_t>,BucketHash> buckets;
        std::vector<std::pair<size_t,size_t>> unique;
        constexpr double join = 2e-6;
        auto node = [&](XYZ position, const Piece &piece) {
            const XY p = xy(position);
            auto feature = mesh_feature(piece,position);
            if (feature[0] == 1) {
                const MeshEdge edge = edge_root({feature[1],feature[2]});
                feature = {1,edge.first,edge.second};
            }
            const int64_t x = std::llround(p.x/join), y = std::llround(p.y/join);
            for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) {
                auto it = buckets.find({feature,x+dx,y+dy});
                if (it != buckets.end()) for (size_t n : it->second)
                    if (std::hypot(nodes[n].p.x-p.x,nodes[n].p.y-p.y) <= join) return n;
            }
            const size_t n = nodes.size();
            nodes.push_back({p,{}}); buckets[{feature,x,y}].push_back(n); return n;
        };
        auto line = [&](XYZ a, XYZ b, const Piece &piece) {
            const size_t na = node(a,piece), nb = node(b,piece);
            if (na != nb) unique.push_back({na,nb});
        };
        for (size_t cell = 0; cell < cells.size(); ++cell) {
            cancel();
            const auto &tessellation = surface.cells[cell];
            if (tessellation.vertices.empty()) continue;
            auto first = std::lower_bound(cells[cell].begin(),cells[cell].end(),tessellation.low+offset,
                [](const Piece &p, double z) { return p.max_z < z; });
            const size_t sx = tessellation.sx, sy = tessellation.sy;
            const auto &vertices = tessellation.vertices;
            auto facet = [&](std::array<XYZ,3> triangle) {
                for (XYZ &p : triangle) p.z += offset;
                const double lo = std::min({triangle[0].z,triangle[1].z,triangle[2].z});
                const double hi = std::max({triangle[0].z,triangle[1].z,triangle[2].z});
                const auto &index = z_index[cell];
                size_t n = size_t(std::lower_bound(first,cells[cell].end(),lo,
                    [](const Piece &p, double z) { return p.max_z < z; })-cells[cell].begin());
                while (n < index.min_z.size()) {
                    if (n%z_block == 0 && index.block_min_z[n/z_block] > hi) { n += z_block; continue; }
                    if (index.min_z[n] <= hi) trace(cells[cell][n],triangle,line);
                    ++n;
                }
            };
            for (size_t v = 0; v < sy; ++v) for (size_t u = 0; u < sx; ++u) {
                const size_t n = v*(sx+1)+u;
                facet({vertices[n],vertices[n+1],vertices[n+sx+2]});
                facet({vertices[n],vertices[n+sx+2],vertices[n+sx+1]});
            }
        }
        // Shared surface facets may report a segment twice. A tangential mesh
        // edge contributes opposite directions, which cancel. Both operations
        // are scoped by mesh topology through the endpoint identities above.
        // Sorted and deduplicated: the same order the former std::set gave.
        std::sort(unique.begin(),unique.end());
        unique.erase(std::unique(unique.begin(),unique.end()),unique.end());
        for (const auto &[a,b] : unique) {
            if (std::binary_search(unique.begin(),unique.end(),std::make_pair(b,a))) continue;
            nodes[a].edges.push_back(edges.size());
            edges.push_back({a,b});
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
                current = edge.b;
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
    static std::array<size_t,3> mesh_feature(const Piece &piece, XYZ p)
    {
        constexpr double epsilon = 1e-9;
        for (size_t n = 0; n < 3; ++n)
            if (distance(p,piece.polygon[n]) <= epsilon)
                return {0,piece.vertex_ids[n],0};
        // Thin triangles can put an endpoint within epsilon of two edges.
        // Use the nearest edge so the short segment between them is preserved.
        double closest = epsilon*epsilon;
        std::array<size_t,3> feature {2,piece.face_id,0};
        for (size_t n = 0; n < 3; ++n) {
            const size_t next = (n+1)%3;
            const XYZ edge = subtract(piece.polygon[next],piece.polygon[n]);
            const XYZ offset = subtract(p,piece.polygon[n]);
            const double length_squared = dot(edge,edge);
            const XYZ perpendicular = cross(edge,offset);
            const double squared_distance = dot(perpendicular,perpendicular)/length_squared;
            if (squared_distance <= closest &&
                dot(offset,edge) >= 0. && dot(offset,edge) <= length_squared) {
                const auto ids = std::minmax(piece.vertex_ids[n],piece.vertex_ids[next]);
                closest = squared_distance;
                feature = {1,ids.first,ids.second};
            }
        }
        return feature;
    }

    template<class Emit> static void trace(const Piece &piece, const std::array<XYZ,3> &triangle, Emit emit)
    {
        const XYZ origin = piece.polygon.front();
        std::array<double,3> residual;
        for (size_t n = 0; n < 3; ++n) residual[n] = dot(subtract(triangle[n],origin),piece.normal);
        if (*std::min_element(residual.begin(),residual.end()) > 0. ||
            *std::max_element(residual.begin(),residual.end()) < 0.) return;
        std::array<XYZ,6> hits;
        size_t hit_count = 0;
        auto add = [&](XYZ p) {
            for (size_t n = 0; n < hit_count; ++n) if (distance(p,hits[n]) < 1e-10) return;
            hits[hit_count++] = p;
        };
        for (size_t n = 0; n < 3; ++n) {
            const size_t next = (n+1)%3;
            if (std::abs(residual[n]) < 1e-12) add(triangle[n]);
            if ((residual[n] >= 0.) != (residual[next] >= 0.))
                add(mix(triangle[n],triangle[next],residual[n]/(residual[n]-residual[next])));
        }
        if (hit_count != 2) return;
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
            emit(a, b, piece);
        }
    }
};

} // namespace

Surface Job::surface_at(double z) const
{
    const LayerStack &stack = this->stack();
    const size_t last = stack.layers.size() - 1;
    const Surface *below, *above;
    double fraction, extension = 0., below_offset = 0.;
    if (z <= 0.) {
        below = above = &stack.bed;
        fraction = 0.; extension = z;
    } else if (z >= layer_z(last)) {
        below = above = &stack.layers[last];
        fraction = 0.; extension = z - layer_z(last);
    } else {
        const size_t layer = size_t(std::max(0., std::ceil((z - first_layer_height) / layer_height)));
        const auto offset = stack.references.empty() ? std::optional<double>() : stack.reference_offsets[layer];
        below = stack.references.empty() ? (layer ? &stack.layers[layer - 1] : &stack.bed) :
                offset ? &stack.layers[layer] : &stack.references[layer];
        below_offset = offset.value_or(0.);
        above = &stack.layers[layer];
        fraction = (z - (layer ? layer_z(layer - 1) : 0.)) / (layer ? layer_height : first_layer_height);
    }
    Surface surface;
    surface.plane = false;
    surface.patches.reserve(below->patches.size());
    for (size_t cell = 0; cell < below->patches.size(); ++cell) {
        Patch patch = below->patches[cell];
        for (double &height : patch.control) height += below_offset;
        for (size_t n = 0; n < patch.control.size(); ++n)
            patch.control[n] += fraction * (above->patches[cell].control[n] - patch.control[n]) + extension;
        surface.patches.push_back(std::move(patch));
    }
    return surface;
}

Vec3d Job::to_model(const Vec3d &point) const
{
    const LayerStack &stack = this->stack();
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
        const auto offset = stack.references.empty() ? std::optional<double>() : stack.reference_offsets[layer];
        const double below = stack.value(stack.references.empty() ? (layer ? stack.layers[layer - 1] : stack.bed) :
                                         offset ? stack.layers[layer] : stack.references[layer], xy) + offset.value_or(0.);
        z = below + fraction * (stack.value(stack.layers[layer], xy) - below);
    }
    return {point.x(), point.y(), z};
}

std::optional<double> Job::surface_offset(double z) const
{
    const LayerStack &stack = this->stack();
    if (z <= 0.) return std::nullopt;
    if (z >= layer_z(stack.layers.size()-1)) return z-layer_z(stack.layers.size()-1);
    const size_t layer = size_t(std::max(0.,std::ceil((z-first_layer_height)/layer_height)));
    if (stack.references.empty() || !stack.reference_offsets[layer]) return std::nullopt;
    const double fraction = (z-(layer ? layer_z(layer-1) : 0.))/(layer ? layer_height : first_layer_height);
    return *stack.reference_offsets[layer]*(1.-fraction);
}

struct MeshSlicer::Impl {
    using Clock = std::chrono::steady_clock;
    const Job &job;
    std::unique_ptr<Intersection> intersection;
    // Meshes may be requested in multiple Z ranges. Cache only final contours,
    // with identical float Z keys for printable and reference calls.
    std::map<float,std::array<Polygons,3>> contours;
    std::mutex mutex;
    double index_time = 0., contour_time = 0.;

    Impl(const Job &j, const Object &object, const std::function<void()> &cancel) : job(j)
    {
        const auto started = Clock::now();
        for (const auto &triangle : object.triangles) for (const auto &p : triangle.vertices)
            if (!job.stack().covers(xy(p))) invalid(_u8L("shared grid does not cover the placed model footprint"));
        intersection = std::make_unique<Intersection>(job.stack(),object,cancel);
        index_time = std::chrono::duration<double>(Clock::now()-started).count();
    }
};

MeshSlicer::MeshSlicer(const Job &job, const Object &object, const std::function<void()> &cancel)
    : impl(std::make_unique<Impl>(job,object,cancel)) {}
MeshSlicer::~MeshSlicer() = default;
double MeshSlicer::index_seconds() const { return impl->index_time; }
double MeshSlicer::contour_seconds() const { return impl->contour_time; }

std::vector<Polygons> MeshSlicer::cross_sections(const std::vector<float> &zs, const std::function<void()> &cancel,
                                               double reference_offset, bool reference)
{
    std::lock_guard<std::mutex> lock(impl->mutex);
    const Job &job = impl->job;
    const bool masked = job.stack().masked;
    std::vector<float> missing;
    for (float z : zs) if (!masked || !impl->contours.count(z)) missing.push_back(z);
    std::vector<std::array<Polygons,3>> computed(missing.size());
    const auto started = Impl::Clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,missing.size(),1), [&](const tbb::blocked_range<size_t> &range) {
        for (size_t n = range.begin(); n < range.end(); ++n) {
            cancel();
            const float z = missing[n];
            const size_t k = std::min(job.stack().layers.size()-1,
                size_t(std::max(0.,std::ceil((double(z)-job.first_layer_height)/job.layer_height))));
            auto slice = [&](const Intersection::Tessellation &surface, double offset, Polygons &polygons) {
                for (const auto &contour : impl->intersection->cross_section(surface,offset,cancel)) {
                    Points points;
                    points.reserve(contour.size());
                    for (const auto &p : contour) points.emplace_back(p.x,p.y);
                    polygons.emplace_back(std::move(points));
                }
            };
            const auto offset = job.surface_offset(z);
            if (masked) {
                // One layer task owns its tessellation and releases it before
                // taking another layer, rather than retaining N grid meshes.
                const auto top = impl->intersection->tessellate(job.stack().layers[k],
                    std::min(-job.layer_height,offset.value_or(0.)),
                    std::max(job.layer_height,offset.value_or(0.)),cancel);
                slice(top,-job.layer_height,computed[n][1]);
                slice(top,job.layer_height,computed[n][2]);
                if (offset) slice(top,*offset,computed[n][0]);
                else slice(impl->intersection->tessellate(job.surface_at(z),0.,0.,cancel),0.,computed[n][0]);
            } else {
                slice(impl->intersection->tessellate(job.surface_at(z),0.,0.,cancel),0.,computed[n][0]);
            }
        }
    });
    impl->contour_time += std::chrono::duration<double>(Impl::Clock::now()-started).count();
    if (!masked) {
        std::vector<Polygons> result;
        result.reserve(computed.size());
        for (auto &layer : computed) result.push_back(std::move(layer[0]));
        return result;
    }
    for (size_t n = 0; n < missing.size(); ++n) impl->contours.emplace(missing[n],std::move(computed[n]));
    std::vector<Polygons> result;
    result.reserve(zs.size());
    const size_t pass = reference ? (reference_offset < 0. ? 1 : 2) : 0;
    for (float z : zs) result.push_back(impl->contours.at(z)[pass]);
    return result;
}
} // namespace Slic3r::Nonplanar
