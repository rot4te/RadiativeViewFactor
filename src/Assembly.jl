# src/Assembly.jl
module Assembly

using LinearAlgebra
using SparseArrays

import ..MeshIO:           MeshData, SurfaceElement
import ..Quadrature:       gauss_legendre_2d
import ..BVH:              BVHTree, build_bvh
import ..ViewFactorKernel: element_pair_view_factor

export ViewFactorResult,
       compute_view_factors,
       aggregate_by_group,
       check_reciprocity,
       check_closure

struct ViewFactorResult
    F_elem      :: Matrix{Float64}
    A_elem      :: Vector{Float64}
    F_group     :: Matrix{Float64}
    A_group     :: Vector{Float64}
    group_tags  :: Vector{Int}
    group_names :: Vector{String}
end

"""
    compute_view_factors(mesh; nquad=4, obstruction_groups=Int[],
                         self_vf=false, verbose=true) -> ViewFactorResult

Compute all element-pair view factors.

Arguments
---------
- `mesh`               : `MeshData` from `load_mesh`
- `nquad`              : Gauss points per direction (nquad² points per element pair)
- `obstruction_groups` : vector of physical group tags that act as potential
                         obstructors. For each emitting/receiving group pair,
                         the tags of those two groups are removed from this
                         vector before building the BVH, so a surface never
                         blocks rays that originate from or terminate on it.
                         Pass an empty vector (default) to skip obstruction checking.
- `self_vf`            : compute self view factors (for curved surfaces)
- `verbose`            : print progress
"""
function compute_view_factors(mesh              ::MeshData;
                               nquad            ::Int       = 4,
                               obstruction_groups::Vector{Int} = Int[],
                               self_vf          ::Bool      = false,
                               verbose          ::Bool      = true)::ViewFactorResult

    elems  = mesh.surface_elems
    coords = mesh.coords
    N      = length(elems)

    check_obs = !isempty(obstruction_groups)

    # Pre-build a BVH for each unique set of obstruction groups that will be
    # needed. The key is the sorted tuple of group tags after removing the
    # source and destination groups. Cache so identical sets share one BVH.
    bvh_cache = Dict{Vector{Int}, Union{BVHTree,Nothing}}()

    function get_bvh(group_i::Int, group_j::Int)::Union{BVHTree,Nothing}
        check_obs || return nothing
        active = sort(filter(g -> g != group_i && g != group_j, obstruction_groups))
        isempty(active) && return nothing
        get!(bvh_cache, active) do
            # Concatenate triangle soups for all active obstruction groups
            soups = [mesh.group_tri_soup[g] for g in active if haskey(mesh.group_tri_soup, g)]
            isempty(soups) && return nothing
            total = sum(size(s, 3) for s in soups)
            merged = Array{Float64,3}(undef, 3, 3, total)
            t = 0
            for s in soups
                nt = size(s, 3)
                merged[:, :, t+1:t+nt] .= s
                t += nt
            end
            build_bvh(merged)
        end
    end

    verbose && println("Computing view factors for $N elements (nquad=$nquad)…")
    check_obs && verbose &&
        println("  Obstruction groups: ", [mesh.group_tags[g] for g in obstruction_groups])

    raw_integral = zeros(Float64, N, N)
    A_elem       = zeros(Float64, N)

    # Area computation (no obstruction needed)
    for i in 1:N
        _, Ai     = element_pair_view_factor(coords, elems[i], elems[i],
                                              nquad, nothing)
        A_elem[i] = Ai
    end

    Threads.@threads for i in 1:N
        gi      = elems[i].group
        j_start = self_vf ? i : i + 1
        for j in j_start:N
            gj   = elems[j].group
            bvh  = get_bvh(gi, gj)
            integ, _ = element_pair_view_factor(coords, elems[i], elems[j],
                                                 nquad, bvh)
            raw_integral[i, j] = integ
            raw_integral[j, i] = integ
        end
        verbose && i % max(1, N÷10) == 0 && println("  … row $i / $N done")
    end

    F_elem = raw_integral ./ reshape(A_elem, N, 1)

    group_tags, group_names, F_group, A_group = _aggregate(mesh, F_elem, A_elem)

    if verbose
        println("Done.")
        println("  Row-sum check (element level) — max |Σⱼ Fᵢⱼ - 1| : ",
                maximum(abs.(vec(sum(F_elem, dims=2)) .- 1.0)))
        println("  Row-sum check (group level)   — max |Σⱼ Fᵢⱼ - 1| : ",
                maximum(abs.(vec(sum(F_group, dims=2)) .- 1.0)))
    end

    return ViewFactorResult(F_elem, A_elem, F_group, A_group,
                             group_tags, group_names)
end

function _aggregate(mesh, F_elem, A_elem)
    gtags  = sort(collect(keys(mesh.group_tags)))
    gnames = [mesh.group_tags[t] for t in gtags]
    G      = length(gtags)

    A_group = zeros(Float64, G)
    for (k, tag) in enumerate(gtags)
        for ei in mesh.group_elems[tag]
            A_group[k] += A_elem[ei]
        end
    end

    F_group = zeros(Float64, G, G)
    for (gi, tagi) in enumerate(gtags)
        for ei in mesh.group_elems[tagi]
            for (gj, tagj) in enumerate(gtags)
                Σ = sum(F_elem[ei, ej] for ej in mesh.group_elems[tagj])
                F_group[gi, gj] += A_elem[ei] * Σ
            end
        end
        F_group[gi, :] ./= A_group[gi]
    end

    return gtags, gnames, F_group, A_group
end

function aggregate_by_group(result::ViewFactorResult, mesh::MeshData)
    tags, names, Fg, Ag = _aggregate(mesh, result.F_elem, result.A_elem)
    return Fg, Ag, tags, names
end

function check_reciprocity(result::ViewFactorResult; tol::Float64=1e-4)::Bool
    F = result.F_elem;  A = result.A_elem;  N = size(F, 1)
    max_err = 0.0
    for i in 1:N, j in i+1:N
        err = abs(A[i]*F[i,j] - A[j]*F[j,i]) / max(A[i]*F[i,j], 1e-30)
        max_err = max(max_err, err)
    end
    println("Reciprocity max relative error: $max_err")
    return max_err < tol
end

function check_closure(result::ViewFactorResult; tol::Float64=1e-3)::Bool
    row_sums = vec(sum(result.F_elem, dims=2))
    println("Row sums: min=$(round(minimum(row_sums),digits=6)), " *
            "max=$(round(maximum(row_sums),digits=6))")
    return maximum(row_sums) <= 1.0 + tol
end

end # module Assembly
