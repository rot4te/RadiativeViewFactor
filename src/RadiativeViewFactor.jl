# src/RadiativeViewFactor
module RadiativeViewFactor

using LinearAlgebra
using StaticArrays
using SparseArrays

include("MeshIO.jl")
include("Quadrature.jl")
include("Geometry.jl")
include("BVH.jl")
include("RayCast.jl")
include("ViewFactorKernel.jl")
include("Assembly.jl")

using .MeshIO:    load_mesh, MeshData
using .Assembly:  compute_view_factors, aggregate_by_group,
                  check_reciprocity, check_closure, ViewFactorResult

export load_mesh,
       compute_view_factors,
       aggregate_by_group,
       check_reciprocity,
       check_closure,
       MeshData,
       ViewFactorResult

end # module
