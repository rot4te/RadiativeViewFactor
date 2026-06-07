# RadiativeViewFactor.jl

[![CI](https://github.com/rot4te/RadiativeViewFactor/actions/workflows/CI.yml/badge.svg)](https://github.com/rot4te/RadiativeViewFactor/actions/workflows/CI.yml)

A Julia package for computing **radiative view factors** between arbitrary surfaces
or curves discretized on 2nd-order meshes generated with [Gmsh](https://gmsh.info/).

## Features

- Reads Gmsh `.msh` files (v2.2 and v4) via the Gmsh Julia SDK
- **3D surface meshes** (`surface_dim=2`): Quad8, Quad9 (centre node dropped), and Tri6 elements
- **2D planar curve meshes** (`surface_dim=1`): Line3 elements; computes view factors per unit depth using the 2D kernel cos θᵢ cos θⱼ / (2r)
- Groups radiating geometry by **Gmsh physical groups**; view factors reported at both element and group level; rows and columns of `F_group` indexed by `result.group_tags` / `result.group_names`
- **Three integration methods** selectable per call:
  - *Gauss–Legendre quadrature* (default): pre-tabulated for n ≤ 5, Golub–Welsch algorithm for n > 5; Dunavant rules for triangular elements; 1-D Gauss–Legendre for Line3 curve elements
  - *Monte Carlo*: stratified area sampling with O(1/N) variance convergence per element pair; per-thread independent RNG streams on CPU; xorshift64 per-thread PRNG on GPU
  - *Duffy transformation* (`use_duffy=true`): Sauter–Schwab singularity regularisation for Quad8 element pairs sharing a vertex (8-region decomposition) or edge (5-region decomposition); falls back to standard quadrature for non-adjacent pairs and non-Quad8 families
- **Obstruction detection** via ray–triangle (3D) or ray–segment (2D) intersection on an axis-aligned BVH; works on both CPU and GPU backends and with all three integration methods
- `obstruction_groups` interface: pass physical group tags of potential occluders; source and destination groups are automatically excluded per pair
- Automatic **normal orientation correction** at load time: for curve meshes (`surface_dim=1`), element normals are oriented to point toward the adjacent transfinite surface interior, determined from mesh element connectivity (works with `.msh` v2.2 and v4)
- Optional **normal reversal**: `reverse_normals=true` flips all normals; `reverse_groups=[...]` flips specific physical groups only
- **Mesh visualisation** via an optional Makie extension: load any Makie backend and call `plot_mesh_normals(mesh)` to inspect element geometry and normal directions
- CPU backend: multi-threaded via `Threads.@threads`
- GPU backends: NVIDIA (`CUDABackend`, Float64) and Apple Silicon (`MetalBackend`, Float32) via KernelAbstractions.jl; stackless BVH traversal eliminates per-thread stack memory pressure
- **Reciprocity** and **closure** (row-sum) checks on the assembled matrix

## Project Layout

```
RadiativeViewFactor.jl/
├── src/
│   ├── RadiativeViewFactor.jl   # Package entry-point and public exports
│   ├── MeshIO.jl                # Gmsh mesh loading; element reading; normal orientation
│   ├── Quadrature.jl            # Gauss–Legendre (1-D and 2-D) and Dunavant rules
│   ├── Geometry.jl              # Shape functions, normals, Jacobians for all element types
│   ├── BVH.jl                   # Axis-aligned BVH; triangle and segment soup support
│   ├── RayCast.jl               # CPU visibility test; dispatches on mesh_dim
│   ├── ViewFactorKernel.jl      # 3D and 2D deterministic kernels; element-pair integrator
│   ├── DuffyKernel.jl           # Sauter–Schwab Duffy transformation for singular pairs
│   ├── MCKernel.jl              # CPU Monte Carlo integrator with stratified sampling
│   ├── Results.jl               # ViewFactorResult, _aggregate, check functions
│   ├── GPUBVH.jl                # Stackless flat BVH for GPU: build + inline traversal
│   ├── GPUKernels.jl            # KernelAbstractions deterministic kernels (Quad8 + Tri6)
│   ├── GPUMCKernels.jl          # KernelAbstractions Monte Carlo kernel (xorshift64 PRNG)
│   ├── Assembly.jl              # CPU assembly; integration dispatch; GPU hook registry
│   └── GPUAssembly.jl           # GPU assembly path; registers GPU hook at load time
├── ext/
│   ├── RadiativeViewFactorCUDAExt.jl    # Registers CUDABackend → CuArray, Float64
│   ├── RadiativeViewFactorMetalExt.jl  # Registers MetalBackend → MtlArray, Float32
│   └── RadiativeViewFactorMakieExt.jl  # plot_mesh_normals (any Makie backend)
├── test/
│   └── runtests.jl
└── Project.toml
```

## Quick Start

### 3D surface mesh — deterministic quadrature

```julia
using RadiativeViewFactor

mesh   = load_mesh("geometry.msh")   # surface_dim=2 is the default
result = compute_view_factors(mesh; nquad=4)

# F_group[i,j] = view factor from group_names[i] to group_names[j]
println(result.group_names)
println(result.F_group)

# Look up a specific pair by name
i = findfirst(==("emitter"),  result.group_names)
j = findfirst(==("receiver"), result.group_names)
println("F(emitter → receiver) = ", result.F_group[i, j])

check_reciprocity(result)   # prints max relative error of Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ
check_closure(result)       # prints row-sum range
```

### With Duffy transformation (near corner/edge singularities)

```julia
# Automatically detects shared vertices and edges between Quad8 elements
# and applies the Sauter–Schwab regularisation only where needed.
result = compute_view_factors(mesh; nquad=6, use_duffy=true)
```

### Monte Carlo

```julia
result = compute_view_factors(mesh; monte_carlo=true, n_samples=50000)

# Reproducible result with a fixed seed
using Random
result = compute_view_factors(mesh; monte_carlo=true, n_samples=50000,
                               rng=MersenneTwister(42))
```

### 2D planar curve mesh (per unit depth)

```julia
# Physical Curve groups required; Line3 elements (Mesh.ElementOrder = 2)
mesh   = load_mesh("planar.msh"; surface_dim=1)
result = compute_view_factors(mesh; nquad=6)
# F_group values are view factors per unit depth

# Normal orientation is corrected automatically toward the adjacent
# transfinite surface interior. Override if needed:
mesh = load_mesh("planar.msh"; surface_dim=1, reverse_normals=true)
mesh = load_mesh("planar.msh"; surface_dim=1, reverse_groups=[2, 5])
```

### Obstruction detection

```julia
# Tags of physical groups that may block rays.
# Source and destination groups are excluded automatically per pair.
result = compute_view_factors(mesh; nquad=4, obstruction_groups=[3, 4])

# Also works with Monte Carlo and Duffy:
result = compute_view_factors(mesh; monte_carlo=true, n_samples=50000,
                               obstruction_groups=[3, 4])
result = compute_view_factors(mesh; nquad=6, use_duffy=true,
                               obstruction_groups=[3, 4])
```

### GPU — NVIDIA CUDA

```julia
using CUDA
result = compute_view_factors(mesh; nquad=4, backend=CUDABackend())
result = compute_view_factors(mesh; monte_carlo=true, n_samples=50000,
                               backend=CUDABackend(), obstruction_groups=[3])
```

### GPU — Apple Metal

```julia
using Metal
# Metal uses Float32 internally; results are promoted to Float64
result = compute_view_factors(mesh; nquad=4, backend=MetalBackend())
```

### Mesh visualisation

```julia
using Plots
using RadiativeViewFactor

mesh = load_mesh("geometry.msh"; surface_dim=1)
fig  = plot_mesh_normals(mesh)

fig = plot_mesh_normals(mesh;
        normal_scale  = 0.05,
        show_nodes    = true,
        show_indices  = true,
        group_colors  = Dict(1=>:red, 2=>:blue))

save("normals.png", fig)   # requires CairoMakie
```

> **Full API documentation** is available at the [package documentation site](https://rot4te.github.io/RadiativeViewFactor.jl).

## Theory

### 3D view factor (surface meshes)

$$F_{ij} = \frac{1}{A_i} \iint_{A_i} \iint_{A_j} \frac{\cos\theta_i \cos\theta_j}{\pi r^2} \, H_{ij} \, dA_j \, dA_i$$

### 2D view factor (curve meshes, per unit depth)

$$F_{ij} = \frac{1}{L_i} \int_{L_i} \int_{L_j} \frac{\cos\theta_i \cos\theta_j}{2r} \, H_{ij} \, dL_j \, dL_i$$

The factor of 2 rather than π in the denominator follows from integrating the 2D radiation intensity over the hemisphere, which gives π/2 rather than π.

### Gauss–Legendre quadrature

The double integral is evaluated at fixed tensor-product Gauss points on each element pair. Pre-tabulated rules for n ≤ 5 use classical nodes and weights (Abramowitz & Stegun §25.4); the Golub–Welsch algorithm is used for n > 5. Dunavant rules are used for triangular elements. Convergence is spectral for smooth integrands but degrades near corner singularities.

### Duffy transformation

For Quad8 element pairs sharing a vertex or edge, the `1/r²` singularity in the kernel is integrable but not efficiently resolved by standard quadrature. The Duffy transformation introduces a radial coordinate ρ measuring distance to the singular point; the Jacobian ρ³ of the 4D transformation cancels the singularity, leaving a bounded integrand on which Gauss quadrature converges rapidly. The implementation uses the Sauter–Schwab decomposition:

- **Common vertex**: 8-region decomposition, each integrated with `nquad⁴` points (total `8 × nquad⁴`)
- **Common edge**: 5-region decomposition, each integrated with `nquad⁴` points (total `5 × nquad⁴`)

Pairs with no shared nodes use standard quadrature (`nquad⁴` points). The singularity type is detected automatically from shared global node indices.

Note: the 2D kernel `1/r` for Line3 elements produces a logarithmic divergence (not `1/r²`) at shared endpoints, which is physically real and not regularisable by the Duffy transformation. `use_duffy` has no effect for `surface_dim=1`.

### Monte Carlo integration

The MC estimator for each element pair draws N stratified sample pairs (xᵢ, xⱼ):

$$\iint K \, dA_j \, dA_i \approx \frac{A_i \cdot A_j}{N} \sum_{k=1}^{N} K(x_i^{(k)}, n_i^{(k)}, x_j^{(k)}, n_j^{(k)}) \cdot H_{ij}^{(k)}$$

Samples are drawn on a ⌊√N⌋ × ⌊√N⌋ stratified grid within the reference element, giving O(1/N) variance convergence for smooth integrands rather than O(1/√N) for plain Monte Carlo. Near corner singularities the variance of the MC estimator diverges (infinite variance for the `1/r²` kernel), making `use_duffy` preferable for those geometries.

On GPU, each thread uses an independent xorshift64 pseudo-random number stream seeded by mixing the global seed with the thread index via the splitmix64 hash.

### Obstruction detection

**CPU**: a BVH is built once per unique set of active obstruction groups and reused across all pairs. Triangle soup for 3D (Möller–Trumbore ray–triangle intersection); line-segment soup for 2D (Cramer's rule line–line intersection).

**GPU**: the BVH is flattened to plain device arrays with miss-link pointers for stackless traversal (no per-thread MVector). Per-triangle group tags allow each GPU thread to skip triangles belonging to the emitter or receiver group without host-side pre-filtering.

## Mesh Requirements

### Surface meshes (`surface_dim=2`)

- 2nd-order elements: `Quad8` (type 16), `Quad9` (type 10, centre node dropped), or `Tri6` (type 9)
- Radiating surfaces in **Physical Surface** groups; obstructors in separate Physical Surface groups
- `Mesh.ElementOrder = 2` before meshing

### Curve meshes (`surface_dim=1`)

- `Line3` elements (Gmsh type 8); `Mesh.ElementOrder = 2` before meshing
- Radiating curves in **Physical Curve** groups
- Normal orientation corrected automatically; use `reverse_normals` or `reverse_groups` if needed

```gmsh
Mesh.ElementOrder = 2;
Physical Curve("emitter")     = {1};
Physical Curve("receiver")    = {2};
Physical Curve("obstruction") = {3};
```

## Performance Notes

### Integration method selection

| Method | Best for | Avoid when |
|---|---|---|
| Quadrature | Smooth geometry, well-separated surfaces | Near corner/edge singularities |
| Duffy | Shared vertices/edges between Quad8 elements | GPU, Monte Carlo, Tri6/Line3 elements |
| Monte Carlo | Many obstructions, near-singular pairs, rough estimates | High-accuracy requirements with few obstructions |

### GPU vs CPU crossover

GPU backends outperform CPU for N ≳ 300–500 elements (quadrature) or N ≳ 200 (Monte Carlo). For small meshes, kernel compilation and host↔device transfer dominate. Metal uses Float32; results are promoted to Float64 before aggregation.

### Duffy transformation cost

The Duffy path evaluates `8 × nquad⁴` (vertex) or `5 × nquad⁴` (edge) kernel evaluations per singular pair, compared to `nquad⁴` for standard quadrature. Since only adjacent element pairs trigger Duffy, the overhead is proportional to the number of shared edges/vertices in the mesh. A lower `nquad` (e.g. 4) with `use_duffy=true` typically outperforms a higher `nquad` (e.g. 16) with standard quadrature for inclined-plate geometries.

## References

The following works informed the numerical methods in this package:

**View factor theory and quadrature:**
- Howell, J. R., Mengüç, M. P., & Siegel, R. (2020). *Thermal Radiation Heat Transfer* (7th ed.). CRC Press. — View factor definitions, reciprocity, crossed-string method, and analytical reference cases.
- Hamilton, D. C., & Morgan, W. R. (1952). *Radiant interchange configuration factors*. NACA Technical Note 2836. — Original tabulation of configuration factor formulae.

**Isoparametric finite elements:**
- Zienkiewicz, O. C., Taylor, R. L., & Zhu, J. Z. (2005). *The Finite Element Method: Its Basis and Fundamentals* (6th ed.). Elsevier. — Quad8 serendipity and Tri6 shape functions, isoparametric mapping, Gauss quadrature.

**Duffy transformation and boundary element singularity treatment:**
- Sauter, S. A., & Schwab, C. (2011). *Boundary Element Methods*. Springer. — Sauter–Schwab common-vertex and common-edge decompositions (§5.3), the definitive reference for the 4D Duffy regularisation used in `DuffyKernel.jl`.
- Duffy, M. G. (1982). Quadrature over a pyramid or cube of integrands with a singularity at a vertex. *SIAM Journal on Numerical Analysis*, 19(6), 1260–1262. — Original Duffy transformation paper.

**Gaussian quadrature:**
- Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss quadrature rules. *Mathematics of Computation*, 23(106), 221–230. — Golub–Welsch algorithm used in `Quadrature.jl` for n > 5.
- Dunavant, D. A. (1985). High degree efficient symmetrical Gaussian quadrature rules for the triangle. *International Journal for Numerical Methods in Engineering*, 21(6), 1129–1148. — Dunavant triangle quadrature rules used for Tri6 elements.

**Ray–triangle intersection:**
- Möller, T., & Trumbore, B. (1997). Fast, minimum storage ray/triangle intersection. *Journal of Graphics Tools*, 2(1), 21–28. — Algorithm used in `BVH.jl` for obstruction testing.

**BVH construction and traversal:**
- Wald, I., et al. (2007). *Ray Tracing Gems* — Stackless BVH traversal via miss-link (skip pointer) encoding, used in `GPUBVH.jl`.

**Monte Carlo integration:**
- Pharr, M., Jakob, W., & Humphreys, G. (2023). *Physically Based Rendering: From Theory to Implementation* (4th ed.). MIT Press. — Stratified sampling, variance reduction, and Monte Carlo estimators for light transport integrals.

**GPU random number generation:**
- Marsaglia, G. (2003). Xorshift RNGs. *Journal of Statistical Software*, 8(14). — xorshift64 PRNG used in `GPUMCKernels.jl`.
- Steele, G. L., Lea, D., & Flood, C. H. (2014). Fast splittable pseudorandom number generators. *ACM SIGPLAN Notices*, 49(10). — splitmix64 hash used to derive per-thread seeds.

## Dependencies

| Package | Role |
|---|---|
| `Gmsh` | Mesh file I/O |
| `KernelAbstractions` | Backend-agnostic GPU kernels |
| `StaticArrays` | Stack-allocated vectors for hot-path geometry |
| `LinearAlgebra`, `SparseArrays`, `Random`, `Statistics` | Standard library |
| `CUDA` *(optional weak dep)* | NVIDIA GPU backend |
| `Metal` *(optional weak dep)* | Apple Silicon GPU backend |
| `Makie` *(optional weak dep)* | Mesh visualisation via `plot_mesh_normals` |
