# src/MeshIO.jl
module MeshIO

using LinearAlgebra
using StaticArrays
import Gmsh: gmsh

export MeshData, SurfaceElement, load_mesh

const ELEM_INFO = Dict{Int, NamedTuple}(
     9 => (n_nodes=6, n_corners=3, family=:tri),
    16 => (n_nodes=8, n_corners=4, family=:quad),
    10 => (n_nodes=9, n_corners=4, family=:quad),
)

struct SurfaceElement
    nodes  :: Vector{Int}
    group  :: Int
    family :: Symbol
end

"""
    MeshData

Fields
------
- `coords`         : (3 × N_nodes) coordinate matrix
- `surface_elems`  : vector of SurfaceElement (all radiating surfaces)
- `group_tags`     : Dict tag → name
- `group_elems`    : Dict tag → element indices (into surface_elems)
- `group_tri_soup` : Dict tag → (3, 3, N_tris) triangle soup for that group's
                     elements only. Used to build per-group BVHs for obstruction.
"""
struct MeshData
    coords         :: Matrix{Float64}
    surface_elems  :: Vector{SurfaceElement}
    group_tags     :: Dict{Int, String}
    group_elems    :: Dict{Int, Vector{Int}}
    group_tri_soup :: Dict{Int, Array{Float64,3}}
end

function load_mesh(filename::AbstractString;
                   surface_dim::Int  = 2,
                   verbose    ::Bool = true)::MeshData
    isfile(filename) || error("Mesh file not found: $filename")
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.open(filename)
    try
        coords, tag2idx = _read_nodes()
        group_tags      = _read_physical_groups(surface_dim)
        surface_elems, group_elems =
            _read_surface_elements(surface_dim, group_tags, tag2idx, verbose)
        group_tri_soup  = _build_group_tri_soups(coords, surface_elems, group_elems)
        if verbose
            counts = Dict{Symbol,Int}(:quad=>0, :tri=>0)
            for e in surface_elems; counts[e.family] += 1; end
            println("Loaded $(length(surface_elems)) surface elements ",
                    "(quad: $(counts[:quad]), tri: $(counts[:tri])) ",
                    "in $(length(group_tags)) physical group(s).")
        end
        return MeshData(coords, surface_elems, group_tags, group_elems,
                        group_tri_soup)
    finally
        gmsh.finalize()
    end
end

function _read_nodes()
    node_tags, coords_flat, _ = gmsh.model.mesh.getNodes()
    N      = length(node_tags)
    coords = Matrix{Float64}(reshape(coords_flat, 3, N))
    tag2idx = Dict{Int,Int}(Int(t) => i for (i,t) in enumerate(node_tags))
    return coords, tag2idx
end

function _read_physical_groups(dim::Int)::Dict{Int,String}
    groups = Dict{Int,String}()
    for (d, tag) in gmsh.model.getPhysicalGroups()
        d == dim || continue
        groups[Int(tag)] = gmsh.model.getPhysicalName(d, tag)
    end
    isempty(groups) && @warn "No physical groups of dimension $dim found in mesh."
    return groups
end

function _read_surface_elements(dim, group_tags, tag2idx, verbose)
    surface_elems = SurfaceElement[]
    group_elems   = Dict{Int,Vector{Int}}(tag => Int[] for tag in keys(group_tags))
    type_counts   = Dict{Int,Int}()

    for (gtag, _) in group_tags
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, gtag)
        for ent in entities
            elem_types, elem_tags, node_tags_per_elem =
                gmsh.model.mesh.getElements(dim, ent)
            for (etype, etags, ntags) in zip(elem_types, elem_tags, node_tags_per_elem)
                itype = Int(etype)
                haskey(ELEM_INFO, itype) || continue
                info      = ELEM_INFO[itype]
                n_nodes   = info.n_nodes
                n_elems   = length(etags)
                ntags_mat = reshape(ntags, n_nodes, n_elems)
                type_counts[itype] = get(type_counts, itype, 0) + n_elems
                for k in 1:n_elems
                    raw = ntags_mat[:, k]
                    if info.family == :quad
                        node_idx = [tag2idx[Int(raw[a])] for a in 1:8]
                        push!(surface_elems, SurfaceElement(node_idx, gtag, :quad))
                    else
                        node_idx = [tag2idx[Int(raw[a])] for a in 1:6]
                        push!(surface_elems, SurfaceElement(node_idx, gtag, :tri))
                    end
                    push!(group_elems[gtag], length(surface_elems))
                end
            end
        end
    end

    if isempty(surface_elems)
        all_types = Set{Int}()
        for (d, ent) in gmsh.model.getEntities(dim)
            etypes, _, _ = gmsh.model.mesh.getElements(d, ent)
            union!(all_types, Int.(etypes))
        end
        error("""
No supported 2nd-order surface elements found in physical groups.
Supported types: Tri6 (9), Quad8 (16), Quad9 (10).
Element types present at dimension $dim: $(sort(collect(all_types)))
Common causes:
  • Mesh is 1st order — re-run with `Mesh.ElementOrder = 2` or `-order 2`
  • Physical groups defined on wrong dimension (expected $dim)
  • Surface entities not included in a Physical Surface
""")
    end

    if verbose
        type_names = Dict(9=>"Tri6", 16=>"Quad8", 10=>"Quad9")
        for (t, n) in sort(collect(type_counts))
            t == 10 && @info "Quad9 (type 10) found — centre node dropped, treated as Quad8."
            println("  Element type $(get(type_names,t,string(t))): $n elements")
        end
    end

    return surface_elems, group_elems
end

"""Build a separate triangle soup for each physical group."""
function _build_group_tri_soups(coords     ::Matrix{Float64},
                                 elems      ::Vector{SurfaceElement},
                                 group_elems::Dict{Int,Vector{Int}})
    soups = Dict{Int, Array{Float64,3}}()
    for (gtag, idxs) in group_elems
        n_tris = sum(elems[i].family == :quad ? 2 : 1 for i in idxs)
        soup   = Array{Float64,3}(undef, 3, 3, n_tris)
        t = 0
        for i in idxs
            el = elems[i]
            c  = el.nodes
            v1 = @view coords[:, c[1]]
            v2 = @view coords[:, c[2]]
            v3 = @view coords[:, c[3]]
            t += 1
            soup[:, 1, t] .= v1
            soup[:, 2, t] .= v2
            soup[:, 3, t] .= v3
            if el.family == :quad
                v4 = @view coords[:, c[4]]
                t += 1
                soup[:, 1, t] .= v1
                soup[:, 2, t] .= v3
                soup[:, 3, t] .= v4
            end
        end
        soups[gtag] = soup
    end
    return soups
end

end # module MeshIO
