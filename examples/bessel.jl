### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 0abe9cd0-7a26-11eb-0d6d-f986ed77c036
using SpecialFunctions, HypergeometricFunctions, StaticArrays

# ╔═╡ c009caa0-7a27-11eb-28ef-bb4011cdcd50
using ClenshawCurtisBessel

# ╔═╡ 854a5a16-7a26-11eb-0c4f-39e135c68926
begin
	ClenshawCurtisBessel.standard_bessel_moment(10, 15/2, 2)
end

# ╔═╡ 386fbd2c-7a2a-11eb-051b-13b127ecfc30
ClenshawCurtisBessel.standard_bessel_moment(1, 5, 1)

# ╔═╡ Cell order:
# ╠═0abe9cd0-7a26-11eb-0d6d-f986ed77c036
# ╠═c009caa0-7a27-11eb-28ef-bb4011cdcd50
# ╠═854a5a16-7a26-11eb-0c4f-39e135c68926
# ╠═386fbd2c-7a2a-11eb-051b-13b127ecfc30
