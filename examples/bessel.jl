### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 0abe9cd0-7a26-11eb-0d6d-f986ed77c036
using SpecialFunctions, HypergeometricFunctions, StaticArrays

# ╔═╡ c009caa0-7a27-11eb-28ef-bb4011cdcd50
using ClenshawCurtisBessel

# ╔═╡ 854a5a16-7a26-11eb-0c4f-39e135c68926


# ╔═╡ 5b6aa4ba-7a35-11eb-2487-179cf98512c0
begin
	N = 10
	a = 50
end

# ╔═╡ 386fbd2c-7a2a-11eb-051b-13b127ecfc30
struct BrandersPiessensProblem{T,N,M} <: ClenshawCurtisBessel.OliverProblem{T,N,M} 
	a::T
	ν::T
end

# ╔═╡ a9a4ba04-7a4f-11eb-2900-15d6325c1e82
# begin
# 	bpp = BrandersPiessensProblem{Float64,8,2}(10, 10)
# 	ClenshawCurtisBessel.OliverP(bpp, 9, 1.0)
# end

# ╔═╡ d9adcfb6-7a4e-11eb-2ada-fbdffca3d2bc


# ╔═╡ Cell order:
# ╠═0abe9cd0-7a26-11eb-0d6d-f986ed77c036
# ╠═c009caa0-7a27-11eb-28ef-bb4011cdcd50
# ╠═854a5a16-7a26-11eb-0c4f-39e135c68926
# ╠═5b6aa4ba-7a35-11eb-2487-179cf98512c0
# ╠═386fbd2c-7a2a-11eb-051b-13b127ecfc30
# ╠═a9a4ba04-7a4f-11eb-2900-15d6325c1e82
# ╠═d9adcfb6-7a4e-11eb-2ada-fbdffca3d2bc
