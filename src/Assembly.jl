# src/Assembly.jl
module Assembly

using LinearAlgebra
using SparseArrays

import ..MeshIO:           MeshData, SurfaceElement
import ..BVH:              BVHTree, build_bvh
import ..ViewFactorKernel: element_pair_view_factor

export ViewFactorResult,
       compute_view_factors,
       aggregate_by_group,
       check_reciprocity,
       check_closure

# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

struct ViewFactorResult
    F_elem      :: Matrix{Float64}
    A_elem      :: Vector{Float64}
    F_group     :: Matrix{Float64}
    A_group     :: Vector{Float64}
    group_tags  :: Vector{Int}
    group_names :: Vector{String}
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

"""
    compute_view_factors(mesh; nquad=4, check_obstruction=true,
                         self_vf=false, verbose=true) -> ViewFactorResult

Compute all element-pair view factors for the mesh.

Arguments
---------
- `mesh`               : `MeshData` from `load_mesh`
- `nquad`              : number of Gauss points per direction (nquad² per element)
- `check_obstruction`  : if `true`, use BVH ray casting to detect blocked pairs
- `self_vf`            : if `true`, compute self view factors (curved elements)
- `verbose`            : print progress
"""
function compute_view_factors(mesh::MeshData;
                               nquad            ::Int  = 4,
                               check_obstruction::Bool = true,
                               self_vf          ::Bool = false,
                               verbose          ::Bool = true)::ViewFactorResult

    elems  = mesh.surface_elems
    coords = mesh.coords
    N      = length(elems)
    bvh = check_obstruction ? build_bvh(mesh.tri_soup) : nothing
    verbose && println("BVH built over $(size(mesh.tri_soup,3)) triangles.")
    verbose && println("Integrating $N × $N element pairs (nquad=$nquad)…")

    # raw_integral[i,j] = ∬_Aᵢ ∬_Aⱼ K dAⱼ dAᵢ   (not yet divided by Aᵢ)
    # A_elem[i]         = area of element i
    raw_integral = zeros(Float64, N, N)
    A_elem       = zeros(Float64, N)

    # Compute element areas independently of the pair loop, so every A_elem[i]
    # is set regardless of loop order.
    gtr = mesh.group_tri_ranges

    for i in 1:N
        _, Ai    = element_pair_view_factor(coords, elems[i], elems[i],
                                             nquad, nothing, gtr;
                                             check_obstruction=false)
        A_elem[i] = Ai
    end

    # Exploit kernel symmetry: ∬_Aᵢ∬_Aⱼ K dAⱼ dAᵢ = ∬_Aⱼ∬_Aᵢ K dAᵢ dAⱼ
    # so we compute only i < j and copy to (j,i).
    # Note: F[i,j] = raw[i,j]/A[i]  and  F[j,i] = raw[j,i]/A[j].
    # Both raw[i,j] and raw[j,i] equal the same double integral, so the copy
    # is numerically exact (not an approximation).
    Threads.@threads for i in 1:N
        j_start = self_vf ? i : i + 1
        for j in j_start:N
            integ, _ = element_pair_view_factor(
                coords, elems[i], elems[j], nquad, bvh, gtr;
                check_obstruction=check_obstruction)

            raw_integral[i, j] = integ
            raw_integral[j, i] = integ   # symmetric kernel
        end
        verbose && i % max(1, N÷10) == 0 &&
            println("  … row $i / $N done")
    end

    # F[i,j] = raw[i,j] / A[i]  — divide each ROW i by A_elem[i].
    # In Julia, `M ./ v` where v is a plain vector divides column-wise (each
    # column j gets divided by v[j]).  To divide row-wise we must reshape v
    # into a column vector so broadcasting aligns along rows.
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

# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

function _aggregate(mesh::MeshData,
                     F_elem::Matrix{Float64},
                     A_elem::Vector{Float64})

    gtags  = sort(collect(keys(mesh.group_tags)))
    gnames = [mesh.group_tags[t] for t in gtags]
    G      = length(gtags)

    A_group = zeros(Float64, G)
    for (k, tag) in enumerate(gtags)
        for ei in mesh.group_elems[tag]
            A_group[k] += A_elem[ei]
        end
    end

    # A_g F_{g→h} = Σᵢ∈g Aᵢ Σⱼ∈h F_{i→j}
    F_group = zeros(Float64, G, G)
    for (gi, tagi) in enumerate(gtags)
        for ei in mesh.group_elems[tagi]
            for (gj, tagj) in enumerate(gtags)
                Σ = 0.0
                for ej in mesh.group_elems[tagj]
                    Σ += F_elem[ei, ej]
                end
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

# ---------------------------------------------------------------------------
# Post-processing helpers
# ---------------------------------------------------------------------------

function check_reciprocity(result::ViewFactorResult; tol::Float64=1e-4)::Bool
    F = result.F_elem
    A = result.A_elem
    N = size(F, 1)
    max_err = 0.0
    for i in 1:N, j in i+1:N
        denom   = max(A[i]*F[i,j], 1e-30)
        err     = abs(A[i]*F[i,j] - A[j]*F[j,i]) / denom
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
