using ClenshawCurtisBessel
using ClenshawCurtisBessel: BrandersPiessensProblem, generate_BC, assembleρ, assembleP
using DoubleFloats
using HypergeometricFunctions
using LinearAlgebra
using ArbNumerics
using TimerOutputs

# reset_timer!(ClenshawCurtisBessel.to)
# for i in 1:100
@time begin
	bpp = BrandersPiessensProblem{Double64,8,2}(10.0, 3 + 1/2)
	indexBC, MBC = generate_BC(bpp, 256)
	ρ = assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ

	ref = big"0.0002069511037367724863632484164263304887654"
	(sol[11] .- ref) ./ ref
end
# end
# show(ClenshawCurtisBessel.to)

##

@time begin
	bpp = BrandersPiessensProblem{Double64,8,2}(10.0, 3)
	indexBC, MBC = generate_BC(bpp, 200)
	ρ = assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ

	ref = big"-0.0001955581634668262545942178152189286110836"
	(sol[11] .- ref) ./ ref
end
