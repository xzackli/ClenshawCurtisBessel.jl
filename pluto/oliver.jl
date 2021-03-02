### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 1d1c0a56-7946-11eb-051f-fb214769d3f0
using PyPlot, BandedMatrices, SparseArrays, LinearAlgebra; PyPlot.svg(true);

# ╔═╡ 20390c4e-7981-11eb-16e4-a5f422480175
md"""We're implementing the method of [Oliver 1968](https://link.springer.com/article/10.1007/BF02166688), a method for solving recurrence relations of the form
```math
P_n(i) y(i) + P_{n-1}(i) y(i+1) + \cdots + P_0(i) y(i+n) = R(i)
```
where the recurrence is unstable in both the forward and backwards solutions. Oliver's method transforms an unstable forward recurrence into a boundary value problem. This can be accomplished by solving a linear system,

```math
\mathbf{P} \vec{y} = \vec{\rho}.
```
We consider the general case in which we specify ``n`` initial values, such that the terms ``Y(a+1), Y(a+2), \ldots, Y(a+n-m)`` are specified, and on the other end, ``Y(b-m+1), \ldots, Y(b)`` are specified for some integer ``m \in [0,n]``, the number of boundary values specified on the end of the recurrence. In this case, the solution is ``\vec{y}``,

```math
\vec{y} = [ Y(a + n - m + 1), \cdots , Y(b - m) ]^T.
```

"""

# ╔═╡ ca5b9a8e-7981-11eb-2dec-dfb58c50ef7d
"""
Construct a boundary value problem from a recurrence relation. The functions for coefficients and the right side dispatch on this abstract type.

`N` and `M` correspond to the total number of specified initial boundary conditions, and the number specified on the end of the recurrence, respectively.
"""
abstract type OliverProblem{T,N,M} end

# ╔═╡ c82988f2-7986-11eb-184c-4dc4685efe09
struct Example61161{T,N,M} <: OliverProblem{T,N,M} end

# ╔═╡ 1a9b22ea-7989-11eb-1123-d923b0f891f1
md"""
Consider an example recurrence relation, which we will call `Example61161`
```math
6 y(s) - 11 y(s+1) + 6 y(s+2) - y(s+3) = 0.
```
with initial conditions ``y(1) = 2``, ``y(2)=5``, and ``y(3) = 15``. This recurrence relation has an explicit solution
```math
y(s) = 1 - 2^{s-1} + 2 \cdot 3^{s-1}.
```
In this case, we actually have constant coefficients so ``P_0(i) = -1``, ``P_1(i) = 6``, ``P_2(i) = -11``, and ``P_3(i) = 6``. Further, ``R(i) = 0``.
"""

# ╔═╡ c4023f50-799d-11eb-18a4-f98202a8f420
y(s) = 1 - 2^(s-1) + 2 * 3^(s-1)

# ╔═╡ c00af4c4-7988-11eb-1dee-4f042a18d53b
"""
	OliverP(OP::OliverProblem, s, i)

Computes ``P_s(i)``, the coefficient ``s`` at solution index ``i``.
"""
function OliverP(OP::Example61161, s, i)
	if s == 0
		return -1
	elseif s == 1
		return 6
	elseif s == 2
		return -11
	elseif s == 3
		return 6
	end
end

# ╔═╡ 67cda80a-798b-11eb-020d-7de13497098c
"""
	OliverR(OP::OliverProblem, i)

Computes ``R(i)``, the right side of the recurrence relation.
"""
OliverR(OP::Example61161, i) = 0

# ╔═╡ c87eb55c-798b-11eb-05db-ab39057c8d20
"""
	assembleP{T}(OP::OliverProblem{M,N}, a, b) where {M, N}

Generic description of the ``\\mathbf{P}`` matrix within the linear system of Oliver's method.
"""
function assembleP(OP::OliverProblem{T,N,M}, a, b) where {T,N,M}
	Psize = b - N - a
	𝐏 = BandedMatrix{T}(undef, (Psize, Psize), (N-M,M))
	for row in 1:Psize
		for col in max(1, row-N+M):min(row + M, Psize)
			Δm = col - row
			𝐏[row, col] = OliverP(OP, M - Δm, a + row)
		end
	end
	return 𝐏
end

# ╔═╡ 20107d58-7981-11eb-1823-b3f867edeb85
assembleP(Example61161{Float64,3,1}(), 0, 10)

# ╔═╡ aaa88a3e-7946-11eb-3301-2720cf0f2301
"""
	assembleρ{T}(OP::OliverProblem{M,N}, a, b) where {M, N}

Generic description of the ``\\vec{\\rho}`` vector within the linear system of Oliver's method.
"""
function assembleρ(OP::OliverProblem{T,N,M}, a, b, YBC) where {T,N,M}
	Psize = b - N - a
	ρ = zeros(T, Psize)

	for row in 1:(N-M)
		Σ = zero(T)
		for s in 1:(N - M - row + 1)
			# every row, increase by 1 the sum UB, YBC arg, P arg
			Σ = Σ + OliverP(OP, N - s + 1, a + row) * YBC[a + s + row - 1]
		end
		ρ[row] = OliverR(OP, a + row) - Σ
	end

	rₛ, rₑ = (b - N - M + 1 - a), (b - N - a)
	for row in rₛ:rₑ
		row_from_end = rₑ - row + 1
		Σ = zero(T)
		for s in 0:(M - row_from_end)
			# every row, decrease by 1 the sum UB, YBC arg, P arg
			Σ = Σ + OliverP(OP, s, b - N - row_from_end + 1) *
				YBC[b - s - row_from_end + 1]
		end
		ρ[row] = OliverR(OP, a+row) - Σ
	end

	return ρ
end

# ╔═╡ f29b3902-7946-11eb-3fcc-0be7bcbf2c6c
md"""Recall the initial conditions ``y(1) = 2``, ``y(2)=5``, and ``y(3) = 15``."""

# ╔═╡ 22514dac-799a-11eb-2940-db93f451d9df
begin
	a, b, n, m = 1, 9, 3, 3
	# YBC = Dict([1,2,3] .=> [2,5,15])
	YBC = Dict( (a+1):(b+1) .=> y.((a+1):(b+1)))
	ex = Example61161{Float64,n,m}()
	P = assembleP(ex, a, b)
	ρ = assembleρ(ex, a, b, YBC)

	# """$(P * y.((a+n-m+1):(b-m))) $(ρ)"""
	sol = P \ ρ
	ref = y.((a+n-m+1):(b-m))
	"""$((sol .- ref) ./ ref)"""
end

# ╔═╡ 3e4fe168-7a07-11eb-0a47-0b31536f8773


# ╔═╡ Cell order:
# ╠═1d1c0a56-7946-11eb-051f-fb214769d3f0
# ╟─20390c4e-7981-11eb-16e4-a5f422480175
# ╠═ca5b9a8e-7981-11eb-2dec-dfb58c50ef7d
# ╠═c82988f2-7986-11eb-184c-4dc4685efe09
# ╟─1a9b22ea-7989-11eb-1123-d923b0f891f1
# ╠═c4023f50-799d-11eb-18a4-f98202a8f420
# ╠═c00af4c4-7988-11eb-1dee-4f042a18d53b
# ╠═67cda80a-798b-11eb-020d-7de13497098c
# ╠═c87eb55c-798b-11eb-05db-ab39057c8d20
# ╠═20107d58-7981-11eb-1823-b3f867edeb85
# ╠═aaa88a3e-7946-11eb-3301-2720cf0f2301
# ╟─f29b3902-7946-11eb-3fcc-0be7bcbf2c6c
# ╠═22514dac-799a-11eb-2940-db93f451d9df
# ╠═3e4fe168-7a07-11eb-0a47-0b31536f8773
