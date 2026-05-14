# src/Geometry.jl
# ---------------------------------------------------------------------------
# Isoparametric geometry for the 8-node serendipity quadrilateral (Quad8).
#
# Reference element: (ξ, η) ∈ [-1,1]²
#
# Node numbering (Gmsh convention, 1-based):
#
#   4---7---3
#   |       |
#   8       6
#   |       |
#   1---5---2
# ---------------------------------------------------------------------------

module Geometry

using LinearAlgebra
using StaticArrays

export quad8_shape,
       quad8_physical_point,
       quad8_normal_and_area_element,
       element_area

@inline function quad8_shape(ξ::Float64, η::Float64)
    # Shape functions
    N1 = 0.25*(1-ξ)*(1-η)*(-ξ-η-1)
    N2 = 0.25*(1+ξ)*(1-η)*( ξ-η-1)
    N3 = 0.25*(1+ξ)*(1+η)*( ξ+η-1)
    N4 = 0.25*(1-ξ)*(1+η)*(-ξ+η-1)
    N5 = 0.5*(1-ξ^2)*(1-η)
    N6 = 0.5*(1+ξ)*(1-η^2)
    N7 = 0.5*(1-ξ^2)*(1+η)
    N8 = 0.5*(1-ξ)*(1-η^2)
    N  = SVector(N1,N2,N3,N4,N5,N6,N7,N8)

    # ∂N/∂ξ  (derived via product rule)
    dN1dξ = 0.25*(1-η)*(2ξ+η)
    dN2dξ = 0.25*(1-η)*(2ξ-η)
    dN3dξ = 0.25*(1+η)*(2ξ+η)
    dN4dξ = 0.25*(1+η)*(2ξ-η)
    dN5dξ = -ξ*(1-η)
    dN6dξ =  0.5*(1-η^2)
    dN7dξ = -ξ*(1+η)
    dN8dξ = -0.5*(1-η^2)
    dNdξ  = SVector(dN1dξ,dN2dξ,dN3dξ,dN4dξ,dN5dξ,dN6dξ,dN7dξ,dN8dξ)

    # ∂N/∂η  (derived via product rule)
    dN1dη = 0.25*(1-ξ)*(ξ+2η)
    dN2dη = 0.25*(1+ξ)*(-ξ+2η)
    dN3dη = 0.25*(1+ξ)*(ξ+2η)
    dN4dη = 0.25*(1-ξ)*(-ξ+2η)
    dN5dη = -0.5*(1-ξ^2)
    dN6dη = -(1+ξ)*η
    dN7dη =  0.5*(1-ξ^2)
    dN8dη = -(1-ξ)*η
    dNdη  = SVector(dN1dη,dN2dη,dN3dη,dN4dη,dN5dη,dN6dη,dN7dη,dN8dη)

    return N, dNdξ, dNdη
end

"""
    quad8_physical_point(coords, nodes, ξ, η) -> SVector{3,Float64}

Map (ξ,η) to physical space. `nodes` may be a Vector{Int} or SVector{8,Int}.
"""
@inline function quad8_physical_point(coords::Matrix{Float64},
                                       nodes,
                                       ξ::Float64, η::Float64)::SVector{3,Float64}
    N, _, _ = quad8_shape(ξ, η)
    x = @SVector zeros(3)
    for a in 1:8
        xa = SVector{3,Float64}(coords[1,nodes[a]], coords[2,nodes[a]], coords[3,nodes[a]])
        x  = x + N[a]*xa
    end
    return x
end

"""
    quad8_normal_and_area_element(coords, nodes, ξ, η) -> (n̂, dA)

Compute the unit normal and area element at (ξ,η). `nodes` may be a
Vector{Int} or SVector{8,Int}.
"""
@inline function quad8_normal_and_area_element(coords::Matrix{Float64},
                                                nodes,
                                                ξ::Float64, η::Float64)
    _, dNdξ, dNdη = quad8_shape(ξ, η)
    dxdξ = @SVector zeros(3)
    dxdη = @SVector zeros(3)
    for a in 1:8
        xa   = SVector{3,Float64}(coords[1,nodes[a]], coords[2,nodes[a]], coords[3,nodes[a]])
        dxdξ = dxdξ + dNdξ[a]*xa
        dxdη = dxdη + dNdη[a]*xa
    end
    c  = cross(dxdξ, dxdη)
    dA = norm(c)
    return c/dA, dA
end

"""
    element_area(coords, nodes; nquad=4) -> Float64

Numerically integrate the area of one Quad8 element.
"""
function element_area(coords::Matrix{Float64}, nodes; nquad::Int=4)::Float64
    # Inline a minimal GL rule to avoid a circular module dependency
    pts1d = [-√(3/5), 0.0, √(3/5)]
    wts1d = [5/9, 8/9, 5/9]
    if nquad <= 2
        pts1d = [-1/√3, 1/√3]; wts1d = [1.0, 1.0]
    elseif nquad >= 4
        # Use Golub–Welsch via the parent Quadrature module
        rule = Main.ViewFactors.Quadrature.gauss_legendre_2d(nquad)
        A = 0.0
        for k in 1:size(rule.points, 2)
            ξ, η = rule.points[1,k], rule.points[2,k]
            _, dA = quad8_normal_and_area_element(coords, nodes, ξ, η)
            A += rule.weights[k] * dA
        end
        return A
    end
    # 3-point rule fallback
    A = 0.0
    for (ξ,wξ) in zip(pts1d,wts1d), (η,wη) in zip(pts1d,wts1d)
        _, dA = quad8_normal_and_area_element(coords, nodes, ξ, η)
        A += wξ * wη * dA
    end
    return A
end

end # module Geometry
