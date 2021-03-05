module ClenshawCurtisBessel

using LinearAlgebra, BandedMatrices
using SparseArrays, StaticArrays, OffsetArrays
import HypergeometricFunctions
import ArbNumerics, DoubleFloats, SpecialFunctions

include("specialfunctions.jl")  # define i.e. bessel J for high precision types
include("oliver.jl")  # transform recurrence into a boundary value problem
include("piessensbranders.jl")  # compute modified moments

export momentM

end
