# src/GPUAssembly.jl
# ---------------------------------------------------------------------------
# GPU dispatch path for compute_view_factors.
# Called from Assembly.jl when a non-CPU backend is passed.
# Obstruction checking is CPU-only; if obstruction_groups is non-empty when
# calling with a GPU backend, a warning is issued and it is ignored.
# ---------------------------------------------------------------------------

module GPUAssembly

using LinearAlgebra
using KernelAbstractions

import ..MeshIO:     MeshData
import ..GPUKernels: build_gpu_arrays, launch_vf_kernel!
import ..Assembly:   ViewFactorResult, _aggregate

export compute_view_factors_gpu

"""
    compute_view_factors_gpu(mesh, nquad, backend, FloatT; verbose) -> ViewFactorResult

GPU implementation of compute_view_factors.

`backend`  — a KernelAbstractions backend, e.g. `CUDABackend()` or `MetalBackend()`.
`FloatT`   — element type: `Float64` for CUDA, `Float32` for Metal.
`ArrayT`   — device array constructor, provided by the backend extension.
"""
function compute_view_factors_gpu(mesh    ::MeshData,
                                   nquad  ::Int,
                                   backend,
                                   FloatT ::Type,
                                   ArrayT ;
                                   verbose::Bool = true)::ViewFactorResult
    N = length(mesh.surface_elems)
    verbose && println("GPU compute_view_factors: $N elements, nquad=$nquad, ",
                       "FloatT=$FloatT, backend=$(typeof(backend))")

    # Flatten mesh data and transfer to device
    verbose && print("  Transferring mesh to device… ")
    ga = build_gpu_arrays(mesh, nquad, ArrayT, FloatT)
    verbose && println("done.")

    # Launch kernels
    verbose && print("  Running GPU kernel… ")
    raw_dev, area_dev = launch_vf_kernel!(ga, backend)
    verbose && println("done.")

    # Copy results back to CPU
    raw_cpu  = Array(raw_dev)
    area_cpu = Array(area_dev)

    # Promote to Float64 for all post-processing (aggregation, reciprocity checks)
    raw_f64  = Float64.(raw_cpu)
    area_f64 = Float64.(area_cpu)

    # Divide each row i by A[i] to get F_elem
    F_elem = raw_f64 ./ reshape(area_f64, N, 1)

    group_tags, group_names, F_group, A_group =
        _aggregate(mesh, F_elem, area_f64)

    if verbose
        println("  Row-sum check (element level) — max |Σⱼ Fᵢⱼ - 1| : ",
                maximum(abs.(vec(sum(F_elem, dims=2)) .- 1.0)))
        println("  Row-sum check (group level)   — max |Σⱼ Fᵢⱼ - 1| : ",
                maximum(abs.(vec(sum(F_group, dims=2)) .- 1.0)))
    end

    return ViewFactorResult(F_elem, area_f64, F_group, A_group,
                             group_tags, group_names)
end

end # module GPUAssembly
