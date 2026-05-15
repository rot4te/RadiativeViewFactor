# ext/RadiativeViewFactorCUDAExt.jl
# Loaded automatically when the user does `using CUDA` before or after
# `using RadiativeViewFactor`.  Registers CUDABackend with Float64 and CuArray.
module RadiativeViewFactorCUDAExt

using CUDA
using KernelAbstractions
import .RadiativeViewFactor.Assembly: _gpu_array_type, _gpu_float_type

# NVIDIA GPUs support Float64 natively (on compute-capable hardware)
_gpu_array_type(::CUDABackend) = CuArray
_gpu_float_type(::CUDABackend) = Float64

end # module
