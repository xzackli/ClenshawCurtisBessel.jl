### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 0abe9cd0-7a26-11eb-0d6d-f986ed77c036
using SpecialFunctions, HypergeometricFunctions, StaticArrays

# ╔═╡ c009caa0-7a27-11eb-28ef-bb4011cdcd50
using ClenshawCurtisBessel

# ╔═╡ 2176cb22-7a63-11eb-0e79-377be1de3061
using DoubleFloats

# ╔═╡ d311f31a-7a5f-11eb-3f9b-c3fd0c2655da
using LinearAlgebra

# ╔═╡ 386fbd2c-7a2a-11eb-051b-13b127ecfc30
struct BrandersPiessensProblem{T,N,M} <: ClenshawCurtisBessel.OliverProblem{T,N,M} 
	a::T
	ν::T
end

# ╔═╡ bc00d5ee-7a5d-11eb-2eb5-5b938e7df938
md"""

Computes the modifed moments from Branders & Piessens 1987.

```math
\begin{aligned}
\frac{a^2}{16} &M_{k}(a, \nu) + \left[ (k+1)^2 - \nu^2 - \frac{a^2}{4} \right] M_{k + 2}(a, \nu) + \left[ 4 \nu^2 - 2(k+4) + 4 \right] M_{k+3}(a,\nu) \\
&- \left[ 2(k+4)^2 - 6 + 6\nu^2 - \frac{3a^2}{8} \right] M_{k+4}(a,\nu) + (4\nu^2 + 2(k+4) + 4)M_{k+5}(a,\nu) \\
&+ \left[ (k+7)^2 - \nu^2 - \frac{a^2}{4} \right]M_{k+6} + \frac{a^2}{16} M_{k+8}(a,\nu) \quad = 0
\end{aligned}
```
"""

# ╔═╡ d9adcfb6-7a4e-11eb-2ada-fbdffca3d2bc

function ClenshawCurtisBessel.OliverP(
		OP::BrandersPiessensProblem{T}, s, k) where T
	a, ν = OP.a, OP.ν
	if s == 8
		return a^2/16
	elseif s == 7
		return zero(T)
	elseif s == 6
		return (k+1)^2 - ν^2 - a^2/4
	elseif s == 5
		return 4ν^2 - 2(k+4) + 4
	elseif s == 4
		return -(2 * (k+4)^2 - 6 + 6ν^2 - 3 * a^2 / 8)
	elseif s == 3
		return 4ν^2 + 2(k+4) + 4
	elseif s == 2
		return (k+7)^2 - ν^2 - a^2/4
	elseif s == 1
		return zero(T)
	elseif s == 0
		return a^2/16
	end
end

# ╔═╡ 29239542-7a56-11eb-152d-b386e0588444
"""
	OliverR(OP::BrandersPiessensProblem{T}, i)

Computes ``R(i)``, the right side of the recurrence relation.
"""
ClenshawCurtisBessel.OliverR(
	OP::BrandersPiessensProblem{T}, i) where T = zero(T)

# ╔═╡ dbf33dbe-7a54-11eb-0410-9b91fb2d22c8
import ClenshawCurtisBessel: M₀, M₁, M₂, M₃, Mₖ_asymptotic

# ╔═╡ 97550070-7a54-11eb-3df8-6ff173e946fc
function generate_BC(OP::BrandersPiessensProblem{T,8,2}, maxorder) where T
	a, ν = OP.a, OP.ν
	MBC = Dict{Int,T}()
	MBC[-3] = M₃(a, ν)
	MBC[-2] = M₂(a, ν)
	MBC[-1] = M₁(a, ν)
	MBC[ 0] = M₀(a, ν)
	MBC[ 1] = M₁(a, ν)
	MBC[ 2] = M₂(a, ν)
	MBC[ 3] = M₃(a, ν)
	
	# Branders & Piessens call the max moment N, we call it maxorder
	ℓ = max(maxorder, Int(ceil(a + 10)))
	print(typeof(a), " ", a, " ", ν, " ", ℓ, "\n")
	MBC[ℓ] = Mₖ_asymptotic(a, ν, ℓ)
	MBC[ℓ+1] = Mₖ_asymptotic(a, ν, ℓ+1)
	return ℓ+1, MBC
end

# ╔═╡ 9740e914-7a54-11eb-1a68-0f136447f214
import ClenshawCurtisBessel: assembleP, assembleρ
# using BandedMatrices

# ╔═╡ 3be2f3f0-7a63-11eb-248f-7727012b9b60
HypergeometricFunctions.logabsgamma(x::Double64) = 
	DoubleFloats.loggamma(abs(x)), sign(gamma(x))

# ╔═╡ 6bbd6c44-7a70-11eb-3859-bd57ba21c45d
# SpecialFunctions.besselj(x::Double64, y::AbstractFloat) =
# 	SpecialFunctions.besselj(Float64(x), Float64(x)) 

# ╔═╡ 972d5cdc-7a54-11eb-0619-85688bb73fcf
@time begin
	bpp = BrandersPiessensProblem{DoubleFloat,8,2}(10.0, 3.5)
	indexBC, MBC = generate_BC(bpp, 200)
	ρ = assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ
end

# ╔═╡ ae66d5ca-7a58-11eb-3001-bf4a310c41d4
begin
	ref = big"0.0002069511037367724863632484164263304887654"
	(sol[11] .- ref) ./ ref
end

# ╔═╡ Cell order:
# ╠═0abe9cd0-7a26-11eb-0d6d-f986ed77c036
# ╠═c009caa0-7a27-11eb-28ef-bb4011cdcd50
# ╠═386fbd2c-7a2a-11eb-051b-13b127ecfc30
# ╟─bc00d5ee-7a5d-11eb-2eb5-5b938e7df938
# ╠═d9adcfb6-7a4e-11eb-2ada-fbdffca3d2bc
# ╠═29239542-7a56-11eb-152d-b386e0588444
# ╠═dbf33dbe-7a54-11eb-0410-9b91fb2d22c8
# ╠═97550070-7a54-11eb-3df8-6ff173e946fc
# ╠═9740e914-7a54-11eb-1a68-0f136447f214
# ╠═2176cb22-7a63-11eb-0e79-377be1de3061
# ╠═d311f31a-7a5f-11eb-3f9b-c3fd0c2655da
# ╠═3be2f3f0-7a63-11eb-248f-7727012b9b60
# ╠═6bbd6c44-7a70-11eb-3859-bd57ba21c45d
# ╠═972d5cdc-7a54-11eb-0619-85688bb73fcf
# ╠═ae66d5ca-7a58-11eb-3001-bf4a310c41d4
