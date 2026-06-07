# Mesh Requirements

## Surface meshes (`surface_dim=2`)

Supported 2nd-order element types:

| Gmsh type | Name | Nodes |
|---|---|---|
| 16 | Quad8 | 8 (serendipity quadrilateral) — preferred |
| 10 | Quad9 | 9 (Lagrange quadrilateral) — centre node silently dropped |
| 9  | Tri6  | 6 (quadratic triangle) |

Requirements:
- All radiating surfaces in **Physical Surface** groups
- Obstruction surfaces in separate **Physical Surface** groups
- `Mesh.ElementOrder = 2` set before meshing (or `-order 2` on the command line)

## Curve meshes (`surface_dim=1`)

Supported 2nd-order element types:

| Gmsh type | Name | Nodes |
|---|---|---|
| 8 | Line3 | 3 (quadratic line) |

Requirements:
- All radiating curves in **Physical Curve** groups
- `Mesh.ElementOrder = 2` set before meshing
- Curves must be wound so that element normals point toward opposing surfaces
  (auto-corrected at load time — see below)

## Normal orientation for curve meshes

For `surface_dim=1`, element normals are computed by rotating the tangent
vector 90° counter-clockwise in the xy-plane. The correct orientation depends
on which way the curve is wound in Gmsh.

RadiativeViewFactor.jl corrects normal orientation automatically at load time
by finding the adjacent transfinite surface for each curve entity (from mesh
connectivity, not CAD topology, so `.msh` v2.2 files are supported) and
flipping any element whose normal points away from that surface interior.

If the auto-correction produces wrong results for specific groups:

```julia
# Flip all normals
mesh = load_mesh("planar.msh"; surface_dim=1, reverse_normals=true)

# Flip only specific physical groups (by tag)
mesh = load_mesh("planar.msh"; surface_dim=1, reverse_groups=[2, 5])
```

`reverse_normals` and `reverse_groups` are applied **after** the auto-correction,
so they compose: setting both reverses all normals then re-reverses the specified
groups, with a net effect of leaving those groups alone and flipping all others.

Use [`plot_mesh_normals`](@ref) to inspect the normal directions visually before
running a computation.

## Example Gmsh script

```gmsh
// 2D curve mesh example
Mesh.ElementOrder = 2;

Physical Curve("emitter")     = {1};
Physical Curve("receiver")    = {2};
Physical Curve("obstruction") = {3};
```

```gmsh
// 3D surface mesh example
Mesh.ElementOrder    = 2;
Mesh.RecombineAll    = 1;   // produces Quad8/Quad9 instead of Tri6
Mesh.Algorithm       = 8;   // Frontal-Delaunay for quads

Physical Surface("hotplate")  = {1};
Physical Surface("coldplate") = {2};
Physical Surface("fin")       = {3};
```
