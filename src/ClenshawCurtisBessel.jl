module ClenshawCurtisBessel

using LinearAlgebra, BandedMatrices
using SparseArrays, StaticArrays
import HypergeometricFunctions
import ArbNumerics, DoubleFloats, SpecialFunctions


include("specialfunctions.jl")
include("oliver.jl")
include("bessel.jl")


end
