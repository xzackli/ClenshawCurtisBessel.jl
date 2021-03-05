using FastGaussQuadrature
using BenchmarkTools
import DoubleFloats: Double64
using HypergeometricFunctions: drummond2F0

function L₀truncated(x::T, μ) where T
    N = 30
	Σ = zero(T)
	for n in 1:N
        s = one(T)
        for i in 0:n-1
            s *= (1/2 + μ + i)/(i+1)
            s *= (1/2 - μ + i) / x
        end
		if iseven(n)
			Σ += s
		else
			Σ -= s
		end
	end
	return one(T) + Σ
end

# U(a, b, z::T) = z^(-a) * pFq([a, one(T) + a - b], T[], -one(T)/z)
function L₀(x::T, μ) where T
    return drummond2F0(1/2+μ, 1/2-μ, -1/x; kmax=10_000Int(ceil(μ)))
end

W₀(μ, z) = exp(-z/2) * L₀(z, μ)


function QS_H1(f, a, b, ω::T, ν, n) where T
    i_ω = im / ω  # this is used a lot
    x_gl, w_gl = gausslaguerre(n);
    Σa = complex(zero(T))
    for k in 1:n
        s = w_gl[k] * f(a + i_ω * x_gl[k]) *
            L₀(-2im * ω * (a + i_ω * x_gl[k]), ν) / √(a + i_ω * x_gl[k])
        Σa += s
    end
    Σb = complex(zero(T))
    for k in 1:n
        s = w_gl[k] * f(b + i_ω * x_gl[k]) *
            L₀(-2im * ω * (b + i_ω * x_gl[k]), ν) / √(b + i_ω * x_gl[k])
        Σb += s
    end
    diff = i_ω * exp(im * ω * a) * Σa - i_ω * exp(im * ω * b) * Σb
    res = √(2 / π / ω) * exp(-im * π/2 * (ν+1/2)) * diff
    return real(res)
end

function QS_H2(f, a, b, ω::T, ν, n) where T
    i_ω = im / ω  # this is used a lot
    x_gl, w_gl = gausslaguerre(n);
    Σa = complex(zero(T))
    for k in 1:n
        s = w_gl[k] * f(a - i_ω * x_gl[k]) *
            L₀(2im * ω * (a - i_ω * x_gl[k]), ν) / √(a - i_ω * x_gl[k])
        Σa += s
    end
    Σb = complex(zero(T))
    for k in 1:n
        s = w_gl[k] * f(b - i_ω * x_gl[k]) *
            L₀(2im * ω * (b - i_ω * x_gl[k]), ν) / √(b - i_ω * x_gl[k])
        Σb += s
    end
    diff = i_ω * exp(-im * ω * b) * Σb - i_ω * exp(-im * ω * a) * Σa
    res = √(2 / π / ω) * exp(im * π/2 * (ν+1/2)) * diff
    return real(res)
end


##
Tn = Double64
f(z) = z^(-2)
         #  a    b             ω        ν     n
QS_H1(f, Tn(1.0), Tn(2.0),  Tn(100.0), Tn(400.0), 1000)

##
# z = 0.5
# μ = 3.0
# κ = 4.0
# W(k, m, z) = exp(-z/2) * z^(m + 1/2) * U(1/2 + m - k, 1 + 2m, z)
