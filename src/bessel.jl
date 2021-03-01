

"""
    bessel_moment_G(λ, ν, a)

Computes the standard Bessel moment ``G_{\\lambda, \\nu}(a)`` defined as
```math
G_{\\lambda, \\nu}(a) = \\int_0^a t^{\\lambda} J_{\\nu}(t)\\,dt.
```
"""
bessel_moment_G(λ, ν, a) = pFq(SA[(1 + λ + ν)/2],SA[(3+λ+ν)/2, 1+ν], -a^2/4) * exp(
        log(a) * (1 + λ + ν) - log(2) * (1 + ν) +
        loggamma((1+λ+ν)/2) - loggamma((3+λ+ν)/2) - loggamma(1+ν))


# modified moment boundary conditions (Piessens & Branders 1983 15-18)
M₀(a, ν) = bessel_moment_G(0, ν, a) / a
M₁(a, ν) = ((2/a) * bessel_moment_G(1, ν, a) - bessel_moment_G(0, ν, a)) / a
M₂(a, ν) = ((8/a^2)*bessel_moment_G(2, ν, a) -
    (8/a)*bessel_moment_G(1, ν, a) + bessel_moment_G(0, ν, a)) / a
M₃(a, ν) = (
    (32/a^3) * bessel_moment_G(3, ν, a) -
    (48/a^2) * bessel_moment_G(2, ν, a) +
    (18/a) * bessel_moment_G(1, ν, a) - bessel_moment_G(0, ν, a)
) / a


function Mₖ_asymptotic(a::T, ν, k) where T

    ϕ0 = SA[
        besselj(ν, a),
        zero(T),
        ((3ν-2) * besselj(ν, a) - 3a * besselj(ν-1, a)) / 2,
        zero(T),
        ( 15a * besselj(ν-1, a) + (4 - 15a^2 + 15ν * (ν-1)) * besselj(ν, a) ) / 4,
        0,
        (21a * (-2 + 5a^2 - 5ν^2) * besselj(ν-1, a) +
            (-8 - 105a^2 * (ν-3) + 21ν * (2 + 5 * (ν-1) * ν)) * besselj(ν, a)) / 8
    ]

    ϕπ = SA[
        -besselj(ν, 0),
        zero(T),
        besselj(ν, 0) - (3a/4) * besselj(ν-1, 0),
        zero(T),
        (-15a^2 * besselj(ν-2, 0) + 60a * besselj(ν-1, 0) +
            2 * (15a^2 - 8) * besselj(ν, 0)) / 16,
        zero(T),
        (-105a^3 * besselj(ν-3, 0) +
            840a^2 * besselj(ν-2, 0) -
            (315a^3 - 1008a) * besselj(ν-1, 0) +
            (64 - 1680a^2) * besselj(ν, 0)
        ) / 64
    ]

    Σ = zero(T)
    for j in 0:3
        Σ = Σ + (-1)^j * ((-1)^k * ϕπ[2j+1] - ϕ0[2j+1]) * k^(-2j-2)
    end
    return Σ / 2
end


@doc raw"""
	ClenshawCurtisBessel.OliverP(
		OP::BrandersPiessensProblem{T}, s, k) where T

Computes the modifed moments from Branders & Piessens 1987.

```math
\begin{aligned}
\frac{a^2}{16} &M_{k}(a, \nu) + \left[ (k+1)^2 - \nu^2 - \frac{a^2}{4} \right] M_{k + 2}(a, \nu) + \left[ 4 \nu^2 - 2(k+4) + 4 \right] M_{k+3}(a,\nu) \\
&- \left[ 2(k+4)^2 - 6 + 6\nu^2 - \frac{3a^2}{8} \right] M_{k+4}(a,\nu) + (4\nu^2 + 2(k+4) + 4)M_{k+5}(a,\nu) \\
&+ \left[ (k+7)^2 - \nu^2 - \frac{a^2}{4} \right]M_{k+6} + \frac{a^2}{16} M_{k+8}(a,\nu) \quad = 0
\end{aligned}
```
"""
function ClenshawCurtisBessel.OliverP(
		OP::BrandersPiessensProblem{T}, s, k) where T
	a, ν = OP.a, OP.ν
	if s == 0
		return a^2/16
	elseif s == 1
		return zero(T)
	elseif s == 2
		return (k+1)^2 - ν^2 - a^2/4
	elseif s == 3
		return 4ν^2 - 2(k+4) + 4
	elseif s == 4
		return -(2 * (k+4)^2 - 6 + 6ν^2 - 3 * a^2 / 8)
	elseif s == 5
		return 4ν^2 + 2(k+4) + 4
	elseif s == 6
		return (k+7)^2 - ν^2 - a^2/4
	elseif s == 7
		return zero(T)
	elseif s == 8
		return a^2/16
	end
end
