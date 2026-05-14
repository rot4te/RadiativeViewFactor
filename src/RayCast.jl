# src/RayCast.jl
module RayCast

using StaticArrays
using LinearAlgebra

import ..BVH: BVHTree, intersect_ray_bvh

export is_visible

"""
    is_visible(bvh, x_i, x_j, skip_start, skip_end) -> Bool

Return `true` if the segment from `x_i` to `x_j` is unobstructed.

Triangles with indices in [skip_start, skip_end] (inclusive, 1-based) are
excluded from the BVH query. Pass the union of the tri_ranges of both the
source and destination physical groups so that neither surface blocks rays
that originate from or terminate on it.
"""
@inline function is_visible(bvh       ::BVHTree,
                             x_i      ::SVector{3,Float64},
                             x_j      ::SVector{3,Float64},
                             skip_start::Int,
                             skip_end  ::Int)::Bool
    d = x_j - x_i
    L = norm(d)
    L < eps() && return true
    hit = intersect_ray_bvh(bvh, x_i, d/L, L;
                             skip_start = skip_start,
                             skip_end   = skip_end)
    return !hit
end

end # module RayCast
