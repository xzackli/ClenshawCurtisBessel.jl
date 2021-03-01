module ClenshawCurtisBessel

using LinearAlgebra, BandedMatrices
using SparseArrays, StaticArrays
using SpecialFunctions, HypergeometricFunctions

include("oliver.jl")
include("bessel.jl")

export OliverProblem, assembleP, assembleρ

end
