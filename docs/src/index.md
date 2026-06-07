# RadiativeViewFactor.jl

A Julia package for computing **radiative view factors** between arbitrary surfaces
or curves discretized on 2nd-order meshes generated with [Gmsh](https://gmsh.info/).

## Overview

RadiativeViewFactor.jl evaluates the double-surface integral that defines the
geometric view factor between pairs of finite surfaces, using the mesh generated
by Gmsh as its geometric description. Three integration strategies are available:

- **Gauss–Legendre quadrature** — the default; spectral convergence for smooth
  geometries
- **Monte Carlo** — stratified area sampling; advantageous when many obstructions
  are present or near-singular pairs exist
- **Duffy transformation** — Sauter–Schwab singularity regularization for Quad8
  element pairs sharing a vertex or edge; gives accurate results for inclined
  surfaces with common edges

All three methods support **obstruction detection** via a BVH-accelerated ray
casting, and all work on CPU and GPU backends (Duffy is CPU-only).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/rot4te/RadiativeViewFactor.jl")
```

For GPU support, install the relevant backend before loading the package:

```julia
Pkg.add("CUDA")   # NVIDIA
Pkg.add("Metal")  # Apple Silicon
```

For mesh visualisation:

```julia
Pkg.add("GLMakie")   # interactive window
Pkg.add("CairoMakie") # file output (PNG, SVG, PDF)
```

## Quick Example

```julia
using RadiativeViewFactor

mesh   = load_mesh("geometry.msh")
result = compute_view_factors(mesh; nquad=4)

i = findfirst(==("emitter"),  result.group_names)
j = findfirst(==("receiver"), result.group_names)
println("F(emitter → receiver) = ", result.F_group[i, j])

check_reciprocity(result)
check_closure(result)
```

## Contents

```@contents
Pages = [
    "manual/getting_started.md",
    "manual/mesh_requirements.md",
    "manual/integration_methods.md",
    "manual/obstruction.md",
    "manual/gpu.md",
    "manual/visualisation.md",
    "manual/performance.md",
    "theory.md",
    "api.md",
    "references.md",
]
Depth = 2
```
