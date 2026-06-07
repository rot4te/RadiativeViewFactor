# Getting Started

## Loading a mesh

RadiativeViewFactor.jl reads Gmsh `.msh` files in both v2.2 and v4 format.
All radiating surfaces must belong to named **Physical Surface** groups (3D) or
**Physical Curve** groups (2D).

```julia
using RadiativeViewFactor

# 3D surface mesh (default)
mesh = load_mesh("geometry.msh")

# 2D planar curve mesh (view factors per unit depth)
mesh = load_mesh("planar.msh"; surface_dim=1)
```

## Computing view factors

```julia
result = compute_view_factors(mesh; nquad=4)
```

This returns a [`ViewFactorResult`](@ref) containing view factors at both the
element level (`F_elem`) and the physical-group level (`F_group`).

## Reading the result

The rows and columns of `F_group` correspond to physical groups in the order
given by `result.group_tags` and `result.group_names`:

```julia
# Print all group names
println(result.group_names)

# Look up a specific pair by name
i = findfirst(==("hotplate"),  result.group_names)
j = findfirst(==("coldplate"), result.group_names)
println("F(hotplate → coldplate) = ", result.F_group[i, j])
```

`F_group[i, j]` is the view factor **from** group `i` **to** group `j`.

## Validation checks

```julia
check_reciprocity(result)   # verifies Aᵢ Fᵢⱼ ≈ Aⱼ Fⱼᵢ
check_closure(result)       # verifies row sums ≤ 1
```

Row sums less than 1 are expected for open geometries where radiation escapes
through open boundaries.
