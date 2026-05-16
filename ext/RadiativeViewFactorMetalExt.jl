# ext/RadiativeViewFactorMetalExt.jl
# Loaded automatically when the user does `using Metal` before or after
# `using RadiativeViewFactor`.  Registers MetalBackend with Float32 and MtlArray.
# Apple GPUs do not support Float64 natively; Float32 is used throughout.
module RadiativeViewFactorMetalExt

using Metal
using KernelAbstractions
import RadiativeViewFactor.Assembly: _gpu_array_type, _gpu_float_type

_gpu_array_type(::MetalBackend) = MtlArray
_gpu_float_type(::MetalBackend) = Float32

end # module
