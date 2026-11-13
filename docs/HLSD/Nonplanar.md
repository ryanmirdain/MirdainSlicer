# Direct curved-layer intersections

Supplied nonplanar layers replace mesh-plane intersection. Orca retains its
ordinary volume and region processing, perimeter/infill generators, seam and
path ordering, travel, retraction, cooling, printer scripts and movement writer.
The original placed mesh is not pre-warped or exported as a distorted mesh.

## Coordinate contract

The historical warp retains X and Y and changes Z. A supplied layer top is the
inverse image of a nominal Orca layer height. Nominal heights use the configured
first-layer height and regular layer height. Between these heights the inverse
interpolates the supplied surface patches. Its XY lookup clamps to the grid;
requested X and Y remain unchanged. Above and below the supplied height range,
Z extends at unit rate, matching the historical height-profile inverse.

`Job::surface_at(z)` evaluates this same inverse for a slicing height. Direct
intersection therefore finds contours of the original triangles at the inverse
image of Orca's slicing plane. `Job::to_model(point)` supplies the corresponding
inverse for exported movements. Both operations use the same surface patches
and nominal heights. This is a height-coordinate warp, not an orthonormal
rotation of each layer plane.

Stack validation requires strict ordering between supplied layer surfaces at grid
nodes and corresponding Hermite patch controls. The first surface may cross the
bed plane; extensions outside printable material need not remain above the bed.

The optional top-level JSON boolean `interpolate` defaults to `true`, retaining
cubic XY surface reconstruction. With `false`, each grid cell is two planar
triangles sharing the lower-left to upper-right diagonal. Supplied heights are
preserved and derivatives are ignored. Slicing and output mapping evaluate the
same triangles. Linear interpolation between nominal layer heights remains active.
Beauty exports `false`; Functional exports `true`.

## Intersection seam

`PrintObject::slice()` creates its usual layers from the supplied schedule and
continues through `slice_volumes()` and normal slice cleanup. The per-volume
intersection dispatch receives the same transformations and Z-range filtering
as ordinary slicing. Nonplanar intersections return contours in the object's
normal centered XY frame. Both intersection implementations then call the same
`slice_mesh_ex` contour-closing, fill-rule and simplification implementation.
Volume Booleans, region assignment, polyhole processing and brim grouping stay
in the ordinary pipeline; geometry is not forced into region zero.

Grid-cell clipping selects candidate source triangles only. Intersection uses
the original triangle planes and boundaries; it does not intersect against the
extra edges created by clipping a triangle into separate grid cells.

Contour assembly retains the input mesh's indexed vertex identities. An
intersection endpoint belongs to a mesh vertex, an indexed mesh edge, or a
source face interior. Coordinate tolerances join endpoints only within that
feature, including joins across the surface tessellation. Separate coincident
mesh boundaries are not welded. Repeated directed segments from shared surface
facets are emitted once; opposite segments at a tangential mesh edge cancel.
Loops follow segment direction to retain the winding supplied by mesh normals.
The ordinary polygon fill processing resolves the resulting closed contours;
assembly does not modify or union the input mesh.

Placed instances are sliced separately because their intersections depend on
position within the supplied field. Planar slice caches cannot store this field.
The existing model metadata preserves the stack and placement provenance; the
stack is not duplicated as a single enormous G-code comment.

## Historical output inverse

Orca exports its complete ordinary G-code first. The nonplanar inverse then
maps model motions from the first feature after the first object-start comment
to the explicit machine-end marker. Startup, initial approach and shutdown stay
unchanged. A Z-only transition enters mapped coordinates at fixed XY.

The inverse accepts G0/G1 and XY-plane G2/G3. Moves are subdivided with the
historical maximum-segment-length rule; arcs become linear segments. Total E
and feed rates are retained. There is no local-gap extrusion compensation,
feed-rate rewriting, custom clearance-hop planner or forced lift type. Mapped
moves use absolute XYZ and relative E; original modes and the source extrusion
position are restored for shutdown. The historical 0.01 mm minimum mapped Z is
retained. The regular G-code processor reads the resulting physical movements
for preview and time estimation.

Surface lookup restores Orca's active extruder XY offset before evaluating the
inverse. Commanded XY stays unchanged: the field is in model/nozzle coordinates,
while exported G-code has the extruder offset subtracted. Filament tool commands
select the corresponding configured extruder offset.

The inverse pass is idempotent: its marker prevents re-mapping an already
exported file when Orca reuses a completed export.
