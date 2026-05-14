# src/BVH.jl
# ---------------------------------------------------------------------------
# A simple, self-contained Axis-Aligned Bounding-Volume Hierarchy (BVH) for
# the triangle soup that represents obstruction geometry.
#
# Build once from MeshData.tri_soup; query repeatedly during view-factor
# integration with ray–AABB + ray–triangle tests.
#
# The BVH is a binary tree stored in a flat array (implicit left/right
# children at 2k and 2k+1).  Leaf nodes hold a small list of triangle
# indices; interior nodes hold only the merged AABB.
# ---------------------------------------------------------------------------

module BVH

using StaticArrays
using LinearAlgebra

export BVHTree, build_bvh, intersect_ray_bvh

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

"""Axis-aligned bounding box."""
struct AABB
    lo :: SVector{3, Float64}
    hi :: SVector{3, Float64}
end

"""A node in the BVH tree (stored in a flat vector)."""
struct BVHNode
    aabb       :: AABB
    left       :: Int          # index into node array; 0 = leaf
    right      :: Int
    tri_start  :: Int          # index into sorted triangle list (leaf only)
    tri_count  :: Int
end

"""
    BVHTree

Flat-array BVH over a triangle soup.

Fields
------
- `nodes`     : flat array of BVHNode
- `tri_idx`   : permutation of triangle indices (sorted during build)
- `tri_soup`  : reference to the (3, 3, N) triangle coordinate array
"""
struct BVHTree
    nodes    :: Vector{BVHNode}
    tri_idx  :: Vector{Int}
    tri_soup :: Array{Float64, 3}   # 3 × 3 × N  (vertex, coord, tri)
end

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

const LEAF_MAX = 8   # max triangles per leaf node

"""
    build_bvh(tri_soup) -> BVHTree

Build a BVH over the triangle soup.  `tri_soup` must be (3, 3, N):
  - dim 1: vertex index (1, 2, 3)
  - dim 2: xyz coordinate (1, 2, 3)
  - dim 3: triangle index
"""
function build_bvh(tri_soup::Array{Float64,3})::BVHTree
    N       = size(tri_soup, 3)
    tri_idx = collect(1:N)
    nodes   = BVHNode[]

    _build_recursive!(nodes, tri_soup, tri_idx, 1, N)

    return BVHTree(nodes, tri_idx, tri_soup)
end

function _triangle_aabb(tri_soup::Array{Float64,3}, tidx::Int)::AABB
    lo = SVector(
        min(tri_soup[1,1,tidx], tri_soup[2,1,tidx], tri_soup[3,1,tidx]),
        min(tri_soup[1,2,tidx], tri_soup[2,2,tidx], tri_soup[3,2,tidx]),
        min(tri_soup[1,3,tidx], tri_soup[2,3,tidx], tri_soup[3,3,tidx]),
    )
    hi = SVector(
        max(tri_soup[1,1,tidx], tri_soup[2,1,tidx], tri_soup[3,1,tidx]),
        max(tri_soup[1,2,tidx], tri_soup[2,2,tidx], tri_soup[3,2,tidx]),
        max(tri_soup[1,3,tidx], tri_soup[2,3,tidx], tri_soup[3,3,tidx]),
    )
    return AABB(lo, hi)
end

function _merge_aabb(a::AABB, b::AABB)::AABB
    AABB(min.(a.lo, b.lo), max.(a.hi, b.hi))
end

function _centroid(tri_soup::Array{Float64,3}, tidx::Int)::SVector{3,Float64}
    SVector(
        (tri_soup[1,1,tidx] + tri_soup[2,1,tidx] + tri_soup[3,1,tidx]) / 3,
        (tri_soup[1,2,tidx] + tri_soup[2,2,tidx] + tri_soup[3,2,tidx]) / 3,
        (tri_soup[1,3,tidx] + tri_soup[2,3,tidx] + tri_soup[3,3,tidx]) / 3,
    )
end

"""Recursively build BVH; appends nodes to `nodes` and returns node index."""
function _build_recursive!(nodes::Vector{BVHNode},
                             tri_soup::Array{Float64,3},
                             tri_idx::Vector{Int},
                             start::Int, stop::Int)::Int
    # Compute merged AABB for this range
    box = _triangle_aabb(tri_soup, tri_idx[start])
    for i in start+1:stop
        box = _merge_aabb(box, _triangle_aabb(tri_soup, tri_idx[i]))
    end

    count = stop - start + 1

    # Leaf node
    if count <= LEAF_MAX
        push!(nodes, BVHNode(box, 0, 0, start, count))
        return length(nodes)
    end

    # Choose split axis: longest axis of AABB
    extent = box.hi - box.lo
    axis   = argmax(extent)   # 1=x, 2=y, 3=z

    # Sort triangles by centroid along chosen axis
    mid = (start + stop) ÷ 2
    sort!(view(tri_idx, start:stop);
          by = i -> _centroid(tri_soup, i)[axis])

    # Reserve a slot for this interior node, fill in children later
    push!(nodes, BVHNode(box, 0, 0, 0, 0))
    node_idx = length(nodes)

    left_child  = _build_recursive!(nodes, tri_soup, tri_idx, start, mid)
    right_child = _build_recursive!(nodes, tri_soup, tri_idx, mid+1, stop)

    nodes[node_idx] = BVHNode(box, left_child, right_child, 0, 0)
    return node_idx
end

# ---------------------------------------------------------------------------
# Ray–AABB intersection (slab method)
# ---------------------------------------------------------------------------

"""Return true if ray (origin `o`, direction `d`) hits `box` before t=`tmax`."""
@inline function _ray_aabb(o::SVector{3,Float64},
                            inv_d::SVector{3,Float64},
                            box::AABB, tmax::Float64)::Bool
    t1x = (box.lo[1] - o[1]) * inv_d[1]
    t2x = (box.hi[1] - o[1]) * inv_d[1]
    t1y = (box.lo[2] - o[2]) * inv_d[2]
    t2y = (box.hi[2] - o[2]) * inv_d[2]
    t1z = (box.lo[3] - o[3]) * inv_d[3]
    t2z = (box.hi[3] - o[3]) * inv_d[3]

    tmin = max(min(t1x,t2x), min(t1y,t2y), min(t1z,t2z), 0.0)
    tmax2 = min(max(t1x,t2x), max(t1y,t2y), max(t1z,t2z), tmax)

    return tmin <= tmax2
end

# ---------------------------------------------------------------------------
# Ray–triangle (Möller–Trumbore) — see RayCast.jl for the standalone version
# ---------------------------------------------------------------------------

const _EPS = 1e-10

"""
Return hit distance t > `t_min` if ray hits triangle, else `Inf`.
`v0,v1,v2` are the three vertices of the triangle.
"""
@inline function _ray_triangle(o::SVector{3,Float64},
                                d::SVector{3,Float64},
                                v0::SVector{3,Float64},
                                v1::SVector{3,Float64},
                                v2::SVector{3,Float64},
                                t_min::Float64)::Float64
    e1 = v1 - v0
    e2 = v2 - v0
    h  = cross(d, e2)
    a  = dot(e1, h)
    abs(a) < _EPS && return Inf
    f = 1.0 / a
    s = o - v0
    u = f * dot(s, h)
    (u < 0.0 || u > 1.0) && return Inf
    q = cross(s, e1)
    v = f * dot(d, q)
    (v < 0.0 || u + v > 1.0) && return Inf
    t = f * dot(e2, q)
    t > t_min ? t : Inf
end

# ---------------------------------------------------------------------------
# BVH traversal
# ---------------------------------------------------------------------------

"""
    intersect_ray_bvh(bvh, origin, direction, t_max;
                      skip_start=0, skip_end=0) -> Bool

Return `true` if any triangle in `bvh` is hit by the ray
`origin + t * direction` for `t ∈ (0, t_max)`.

Triangles with 1-based indices in [skip_start, skip_end] are excluded —
pass the group_tri_ranges of both the source and destination groups to
prevent rays from being blocked by the surfaces they originate from or
terminate on.  When the two groups are different their ranges are unioned
into the single contiguous interval [min, max], which may over-exclude a
small number of triangles from third groups if those groups happen to fall
between them in the soup — in practice this is harmless since two surfaces
that can see each other will never have a third surface indexed between them
unless it is genuinely between them geometrically.
"""
function intersect_ray_bvh(bvh::BVHTree,
                             origin::SVector{3,Float64},
                             direction::SVector{3,Float64},
                             t_max::Float64;
                             skip_start::Int = 0,
                             skip_end  ::Int = 0)::Bool
    inv_d = SVector(1.0/direction[1], 1.0/direction[2], 1.0/direction[3])

    stack    = zeros(Int, 64)
    stack[1] = 1
    sp       = 1

    @inbounds while sp > 0
        nidx = stack[sp];  sp -= 1
        node = bvh.nodes[nidx]

        _ray_aabb(origin, inv_d, node.aabb, t_max) || continue

        if node.left == 0   # leaf
            for k in node.tri_start : node.tri_start + node.tri_count - 1
                tidx = bvh.tri_idx[k]
                skip_start <= tidx <= skip_end && continue
                v0 = SVector{3,Float64}(bvh.tri_soup[1,1,tidx],
                                         bvh.tri_soup[1,2,tidx],
                                         bvh.tri_soup[1,3,tidx])
                v1 = SVector{3,Float64}(bvh.tri_soup[2,1,tidx],
                                         bvh.tri_soup[2,2,tidx],
                                         bvh.tri_soup[2,3,tidx])
                v2 = SVector{3,Float64}(bvh.tri_soup[3,1,tidx],
                                         bvh.tri_soup[3,2,tidx],
                                         bvh.tri_soup[3,3,tidx])
                t = _ray_triangle(origin, direction, v0, v1, v2, 0.0)
                t < t_max && return true
            end
        else
            sp += 1;  stack[sp] = node.left
            sp += 1;  stack[sp] = node.right
        end
    end
    return false
end

end # module BVH
