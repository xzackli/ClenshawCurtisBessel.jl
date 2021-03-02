module ClenshawCurtisBessel

using LinearAlgebra, BandedMatrices
using SparseArrays, StaticArrays, OffsetArrays
import HypergeometricFunctions
import ArbNumerics, DoubleFloats, SpecialFunctions

include("specialfunctions.jl")
include("oliver.jl")
include("bessel.jl")

export momentM

end
