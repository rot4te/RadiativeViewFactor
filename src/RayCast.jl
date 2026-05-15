# src/RayCast.jl
module RayCast

using StaticArrays
using LinearAlgebra

import ..BVH: BVHTree, intersect_ray_bvh

export is_visible

"""
    is_visible(bvh, x_i, x_j) -> Bool

Return `true` if the segment from `x_i` to `x_j` is unobstructed by any
triangle in `bvh`.

The BVH should be built from only the triangles that are actual potential
obstructors — i.e. excluding the physical groups of both the emitting and
receiving surfaces. This is handled in Assembly before the BVH is built,
so no skip logic is needed here.
"""
@inline function is_visible(bvh::BVHTree,
                             x_i::SVector{3,Float64},
                             x_j::SVector{3,Float64})::Bool
    d = x_j - x_i
    L = norm(d)
    L < eps() && return true
    return !intersect_ray_bvh(bvh, x_i, d/L, L)
end

end # module RayCast
