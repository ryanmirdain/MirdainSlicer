# Supplied nonplanar layers

The optional native pipeline consumes `orca_nonplanar_layers` version 1 JSON
and the original placed model. Ordinary jobs continue through Orca's standard
PrintObject, region, perimeter, infill and GCode pipeline. A sidecar explicitly
imports a nonplanar job; reopening its saved project preserves that selection.
No model warp, deformed mesh, external service or G-code inverse is involved.

## Coordinates and surfaces

The shared grid is in millimeters in the selected plate's local world XYZ
coordinates, after source units, volume transforms and instance transforms.
Positive Z points upward. `plate_index` is zero-based. The CLI's `--slice` uses
one-based plate numbers; zero means all plates, so it is accepted with a sidecar
only for a single-plate project. The grid must cover every original model vertex
and the entire printable path domain. No evaluation extends beyond its edges.

`plane` means `z = a*x + b*y + c`. `height_grid.z[j][i]` samples the shared
`x[i], y[j]`. Each grid cell is the exact bilinear interpolant, including on
nonuniform grids. Triangulating a cell does not define its surface. Mixed plane
and grid stacks are allowed. IDs are consecutive from zero. H[k] is the nominal
layer TOP/nozzle reference; H[-1] is the supplied bed. The material band is the
original solid intersected with H[k-1] <= z <= H[k].

All coordinates and coefficients must be finite. Axes strictly increase. Runtime
validation checks matrix dimensions, IDs, gaps and interpolated slopes, which
JSON Schema alone cannot express. Gap differences are bilinear, so their extrema
are at cell corners. The gradient norm is convex over each cell and is checked
at all four corners, on both sides of grid boundaries. Vertical gap is distinct
from normal spacing. Normal spacing follows the downward unit normal of H[k]
to its first intersection with H[k-1], solving the piecewise quadratic ray
intersection across grid seams. For parallel local planes this reduces to
`gap / sqrt(1 + |gradient H[k]|^2)`; that formula is not used for general stacks.
The geometry implementation uses 1e-6 mm polygon coordinates and requires
geometry/chord tolerances >= 1e-5 mm. Native export always uses four decimal XYZ
coordinates and requires both tolerances >= 0.0002 mm. Subdivision reserves a
rounding allowance for chord error and maximum segment length. Smaller
tolerances fail explicitly instead of being clamped.

## Original mesh intersection and material regions

`Nonplanar/LayerStack` owns the surface interface (value, gradient, bounds and
cell polynomials). `Nonplanar/Slicer` clips original triangles at shared cell
boundaries. Within a cell, triangle-plane minus layer height is quadratic.
Degree-two Bernstein control values bound it over the complete triangle;
same-sign vertex values do not discard interior crossings or edge tangencies.
Only ambiguous patches refine, by longest-edge bisection. Affine patches clip
directly. Curved boundaries refine to half the geometry tolerance, with
conservative residual margins and a finite recursion/work budget. Exhausting
the budget is an input error, never a planar fallback.

Signed projected face regions use the original oriented mesh's upward ray
winding. Clipper2 positive-winding booleans close contours and preserve islands
and holes. Vertical faces have zero projected area. Coplanar upper faces count
inside; coplanar lower faces are handled by the opposite one-sided limit.
The implementation computes material on BOTH band endpoints and removes the
projection of every boundary strictly inside the band. This catches cavities
that endpoint-only sampling misses. Wall planning uses the full-column part of
that intersection: it conservatively omits partial-height edge columns instead
of extruding a bead outside the original solid.

## Fixed-nozzle bead and flow convention

The nominal line width is a projected XY width. It is intentionally not claimed
to be a constant physical width on a sloped sheet. Walls are round-join inward
offsets of the complete band region by half that width plus geometry tolerance;
subsequent walls step by that projected width. The modeled bead has a vertical
band cross-section with its top attached to H[k] and bottom to H[k-1]. Its full
projected footprint must fit inside the material region. This is a geometric
model, not a guarantee about molten filament or a perfectly smooth finish.

Normal thickness and transverse physical width are stored separately. On a
local patch, transverse width is `projected_width * sqrt(1 + transverse_slope^2)`.
Flow integrates the actual bilinear gap over the full swept rectangular strip,
using exact polygon area/first/mixed moments in each crossed grid cell. Previously
covered footprints on the same band are subtracted before integration, so shared
corners, wall intersections and skin overlaps receive volume once. Flat-ended
infill strips overlap the innermost wall to avoid gaps at their endpoints. Flow is
derived from the newly deposited volume and the final 3D path length. It is
never copied from a planar path. Grid seams split paths before adaptive
segmentation. Chord error, 3D segment length and relative thickness variation
all control segmentation. No automatic normalization, roof pinning or layer
redesign occurs. A final plane matching a flat roof remains that plane.
The preview uses vertical parallelogram sections with the actual transverse
slope. Top and bottom approximation errors are bounded over the entire swept
strip. A line width too wide to represent transverse curvature within tolerance
is explicitly unsupported; reducing segment length never silently flattens it.
Width/height coloring reports physical transverse width and normal spacing,
while rendering uses projected width and vertical gap.

Solid/top/bottom regions follow material presence on neighboring supplied bands
for the configured shell layer counts. Rectilinear infill follows the configured
infill direction and alternates by 90 degrees between layers. Aligned seams use
the lexicographically leftmost contour point, consistently across layers.
Unsupported geometric overhangs and floating islands fail explicitly.

## Machine output, projects and preview

`Nonplanar/Print` adapts placed Print instances and ordinary printer/process/
filament configuration to the native engine. GCodeWriter emits final XYZ/E and
handles retraction and extrusion modes. Every inter-path travel lifts above the
global printed height envelope, bounded across the complete bead width rather
than just its centerline, before XY motion and descent. Fully covered zero-volume
segments use the same clearance travel instead of dragging across the surface.
Each object's first printable band must reach the bed and exceed the extrusion
coordinate resolution. Configured vertical
XYZ/E axis speeds and filament volumetric flow cap deposition speed. Clearance uses
the minimum configured nozzle/rod/lid height conservatively against all already printed
geometry. Failure to find safe Z travel is an error.
Profile feature speeds and first-layer speeds remain in use. A native
deposition-time lower bound drives layer cooling slowdown and the configured
fan thresholds, first-layer fan suppression and fan ramp. Travel adds cooling
time beyond that lower bound. Normal spacing is checked against nozzle profile
height limits. Advanced bridge/support cooling is outside this mode's scope.

Model metadata carries `orca_nonplanar_layers_v1` (the canonical sidecar),
`orca_nonplanar_schema_version`, and placement provenance. The existing 3MF
metadata serializer preserves them in project and G-code 3MF exports. Imported
placement is frozen: changing object/volume/instance transforms requires a new
layer import. G-code carries the same source and schema version, and native
per-segment band dimensions. The processor and libvgcode propagate these through
preview; ordinary moves retain their original representation.

## Input and UI

Schema and examples are in `resources/nonplanar/`. The examples describe a
20 by 20 mm domain and a 1 mm final roof; they are not complete stacks for other
models. Supply layers covering the actual model, rather than scaling examples.

```sh
orca-slicer --load-settings 'printer.json;process.json' \
  --load-filaments filament.json \
  --filament-map-mode Manual --filament-map 1 \
  --arrange 0 --orient 0 --slice 0 \
  --nonplanar-layers layers.json \
  --outputdir out --export-3mf sliced.3mf model.3mf
```

The identical command without the sidecar uses standard slicing for a standard
input project. A saved nonplanar project intentionally retains its imported mode.
GUI: place the original models, select the plate, then File > Import > Import
nonplanar layers. Slice and inspect Preview normally. File > Import > Use standard
layers removes the imported mode. Save Project embeds the input so it can reopen
without the original JSON path. Arrangement, orientation, scaling, duplication,
unit conversion and grounding after import are unsupported.

The native scope uses one tool/filament and one watertight model volume per
object, multiple objects/islands/cavities, classic constant projected-width walls,
rectilinear sparse/solid infill and supplied-layer skins. Set one consistent
`line_width` for every role (or zero for automatic role widths), disable supports,
raft, skirt, brim, prime tower, spiral mode, ironing, fuzzy skin, arc fitting and
overhang speed processing. Set `sparse_infill_pattern` to `rectilinear` and
`print_sequence` to `by layer`. Per-object/per-volume process overrides, modifier
volumes, multiple tools/materials, belt/delta kinematics, slicing plugins, planar
slice caches and G-code postprocessing are unsupported. Unsupported settings
raise specific errors only in nonplanar mode.

Also select `wall_generator=classic` and rectilinear top, bottom and internal
solid patterns. Disable infill combination, machine exclusion regions and
per-layer/timelapse custom scripts. The printable area must be convex so straight
clearance travels between accepted paths stay within its XY limits.
`z_offset` must be zero. Start/end and
filament start/end templates use Orca's placeholder parser; unsupported dynamic
placeholders raise its normal errors instead of being ignored. Templates must
leave a homed fixed-nozzle Cartesian machine; generated moves explicitly restore
absolute XYZ and the configured E mode. The generator waits for the profile's
initial bed/nozzle temperatures and switches to normal temperatures on layer 1.
Select aligned seams without scarf joints or staggered inner seams, and 100%
skin densities. XY/hole/elephant-foot compensation must be zero; filament shrink
compensation must be 100%. These settings would otherwise move original model
material relative to the externally supplied surfaces.
Set `top_shell_thickness` and `bottom_shell_thickness` to zero: this mode uses
shell layer counts and rejects planar minimum-thickness overrides.
