### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 19602f12-7fc4-11eb-33c6-c7fc690e6da2
begin
	using ApproxFun, SpecialFunctions, LinearAlgebra, BenchmarkTools
	import ClassicalOrthogonalPolynomials
end

# ╔═╡ 3497f8a0-7fc4-11eb-09a4-99ede0577168
begin
	using PyPlot
	PyPlot.svg(true)
	function with_pyplot(f::Function)
		f()
		fig = gcf()
		close(fig)
		return fig
	end
end

# ╔═╡ 77a708ca-7fd3-11eb-2c99-3919a290433f
using TimerOutputs

# ╔═╡ 20713f94-7fc4-11eb-051e-fd71664ec440
Tcheb(n,x) = ClassicalOrthogonalPolynomials.chebyshevt(n, x)

# ╔═╡ 2902c6b4-7fc4-11eb-3e3d-171afb2a489c
begin
	S = Chebyshev()
	G = Fun(x -> Tcheb(128, x), S)
	I = ApproxFun.I
end

# ╔═╡ 74d93d2a-7fc4-11eb-0767-f5b765f1f70e
begin
	α = 0.
	ν = 40.
	c = 40
	c2 = Fun(t -> (1+t)^2, S)
	c1 = Fun(t -> (2α + 3)*(1+t), S)
	c0 = Fun(t -> c^2 * (1+t)^2 + (α + 1)^2 - ν^2, S)
	b1 = Fun(t -> (1+t), S)
	b0 = Fun(t -> 2α + 1, S)
	
	bf = BasisFunctional(c/(ν+1) >= 20 ? Integer(ceil(c)) : 0, space(G))
end;

# ╔═╡ 5f850bf2-7fc4-11eb-0b68-437c354291e2
A = [
	c2 * 𝒟^2 + c1 * 𝒟 + c0 * I    0.;
	b1 * 𝒟 + b0 * I               I;       
	bf                            0;
	0                             bf   
];

# ╔═╡ 5f6fa354-7fc4-11eb-3206-dbaa9c57c1fc
b = [-G;
	  0.;
	  0.;
	  0.]

# ╔═╡ 9933d9f4-7fc5-11eb-0472-975929d365b8
F, H = \(A, b, tolerance=1maximum(G.coefficients)*eps())

# ╔═╡ ee258fc0-7fc5-11eb-16fe-7512390173e1
integral(t) = (1+t)^(α+1) * (
	c * (1+t) * (1/2) * (
		besselj(ν-1, c*(1+t)) - besselj(1+ν, c*(1+t))
		) * F(t) + besselj(ν, c*(1+t)) * (α * F(t) + H(t)))

# ╔═╡ a1e75c1c-7fc6-11eb-0b03-ffb83c581358
result = integral(1) - integral(-1)

# ╔═╡ 47101788-803c-11eb-37cc-ef405ac24eed


# ╔═╡ 60a6f6dc-7fc7-11eb-39db-4d3f30fc79aa
begin
	ref = big"-5.331844823289785705118221686111276809121403304e-7"
	"abserr", (result - ref)
end

# ╔═╡ 7ec6de8e-7fc7-11eb-0e82-ab2d6d5cbef6
# @btime \($A, $b, tolerance=10maximum($G.coefficients)*eps())

# ╔═╡ 88cdb6d0-7fd3-11eb-2508-69fe5d96f67a
const to = TimerOutput()

# ╔═╡ fbaa9bac-7fd1-11eb-18a4-d36066740213
function compute_moments(ν, c, α=0.0)
	
	
	S = Chebyshev()
	I = ApproxFun.I
	
	c2 = Fun(t -> (1+t)^2, S)
	c1 = Fun(t -> (2α + 3)*(1+t), S)
	c0 = Fun(t -> c^2 * (1+t)^2 + (α + 1)^2 - ν^2, S)
	b1 = Fun(t -> (1+t), S)
	b0 = Fun(t -> 2α + 1, S)
	
	bf = BasisFunctional(c/(ν+1) >= 20 ? Integer(ceil(c)) : 0, S)
	
	A = [
		c2 * 𝒟^2 + c1 * 𝒟 + c0 * I    0.;
		b1 * 𝒟 + b0 * I               I;       
		bf                            0.;
		0.                             bf   
	];
	
	vals = Float64[]
	qrA = qr(A)  # solve A once
	for n in 0:128
		G = Fun(x -> Tcheb(n, x), S)
		b = [-G;
		  0.;
		  0.;
		  0.
		]
		F, H = \(qrA, b, tolerance=10maximum(G.coefficients)*eps())
		
		t = 1.0
		moment = (1+t)^(α+1) * (
			c * (1+t) * (1/2) * (
			besselj(ν-1, c*(1+t)) - besselj(1+ν, c*(1+t))) * 
			F(t) + besselj(ν, c*(1+t)) * (α * F(t) + H(t)))
		push!(vals, moment)
	end
	return vals
end

# ╔═╡ 69731f42-7fd2-11eb-0c1f-575407f9ee64
@time compute_moments(100.0, 40)

# ╔═╡ b75114de-7fd3-11eb-224e-71b3c667824e
begin
	reset_timer!(to)
	compute_moments(500.0, 300.0)
	to
end

# ╔═╡ b739d9d6-7fd3-11eb-37cc-e7ac8c81c6d9


# ╔═╡ b71f9fa8-7fd3-11eb-3049-05b0c265a291


# ╔═╡ 2700ed36-7fcf-11eb-19c7-59fe0b88eea8
# function do_keller(ν, c, α=0.0)
	
	
# 	S = Chebyshev()
# 	G = Fun(x -> chebyshevt(4, x), S)
# 	I = ApproxFun.I
	
# 	c2 = Fun(t -> (1+t)^2, S)
# 	c1 = Fun(t -> (2α + 3)*(1+t), S)
# 	c0 = Fun(t -> c^2 * (1+t)^2 + (α + 1)^2 - ν^2, S)
# 	b1 = Fun(t -> (1+t), S)
# 	b0 = Fun(t -> 2α + 1, S)
	
# 	bf = BasisFunctional(c/(ν+1) >= 20 ? Integer(ceil(c)) : 0, space(G))
	
# 	A = [
# 		c2 * 𝒟^2 + c1 * 𝒟 + c0 * I    0.;
# 		b1 * 𝒟 + b0 * I               I;       
# 		bf                            0;
# 		0                             bf   
# 	];
	
# 	b = [-G;
# 	  0.;
# 	  0.;
# 	  0.
# 	]
# 	@time F, H = \(A, b, tolerance=10maximum(G.coefficients)*eps())
# 	t = 1.0
	
# 	return (1+t)^(α+1) * (
# 	c * (1+t) * (1/2) * (
# 		besselj(ν-1, c*(1+t)) - besselj(1+ν, c*(1+t))
# 		) * F(t) + besselj(ν, c*(1+t)) * (α * F(t) + H(t)))
# end

# ╔═╡ 8ce980b8-7fcf-11eb-1fbb-fbd0693d81d6
# @time do_keller(ν, c, α, S, G, I, c2, c1, c0, b1, b0, bf, A)

# ╔═╡ 32f414aa-7fd0-11eb-0aa0-efbd6ae3e8af
 # do_keller(3000.0, 1e-1)

# ╔═╡ 65e8471e-7fd0-11eb-0c7a-1550ff5a2094


# ╔═╡ b453481e-7fc5-11eb-2f5e-bb09402f0b88
# with_pyplot() do
# 	# axvline(r)
# 	plot(H.coefficients, "-")
# 	yscale("symlog", linthresh=1e-16)
# end

# ╔═╡ 36cf695a-7fc4-11eb-2968-61eb0cc280bd
# with_pyplot() do
# 	figure(figsize=(10,4))
# 	plot()
# end

# ╔═╡ Cell order:
# ╠═19602f12-7fc4-11eb-33c6-c7fc690e6da2
# ╟─3497f8a0-7fc4-11eb-09a4-99ede0577168
# ╠═20713f94-7fc4-11eb-051e-fd71664ec440
# ╠═2902c6b4-7fc4-11eb-3e3d-171afb2a489c
# ╠═74d93d2a-7fc4-11eb-0767-f5b765f1f70e
# ╠═5f850bf2-7fc4-11eb-0b68-437c354291e2
# ╠═5f6fa354-7fc4-11eb-3206-dbaa9c57c1fc
# ╠═9933d9f4-7fc5-11eb-0472-975929d365b8
# ╠═ee258fc0-7fc5-11eb-16fe-7512390173e1
# ╠═a1e75c1c-7fc6-11eb-0b03-ffb83c581358
# ╠═47101788-803c-11eb-37cc-ef405ac24eed
# ╠═60a6f6dc-7fc7-11eb-39db-4d3f30fc79aa
# ╠═7ec6de8e-7fc7-11eb-0e82-ab2d6d5cbef6
# ╠═77a708ca-7fd3-11eb-2c99-3919a290433f
# ╠═88cdb6d0-7fd3-11eb-2508-69fe5d96f67a
# ╠═fbaa9bac-7fd1-11eb-18a4-d36066740213
# ╠═69731f42-7fd2-11eb-0c1f-575407f9ee64
# ╠═b75114de-7fd3-11eb-224e-71b3c667824e
# ╠═b739d9d6-7fd3-11eb-37cc-e7ac8c81c6d9
# ╠═b71f9fa8-7fd3-11eb-3049-05b0c265a291
# ╠═2700ed36-7fcf-11eb-19c7-59fe0b88eea8
# ╠═8ce980b8-7fcf-11eb-1fbb-fbd0693d81d6
# ╠═32f414aa-7fd0-11eb-0aa0-efbd6ae3e8af
# ╠═65e8471e-7fd0-11eb-0c7a-1550ff5a2094
# ╠═b453481e-7fc5-11eb-2f5e-bb09402f0b88
# ╠═36cf695a-7fc4-11eb-2968-61eb0cc280bd
