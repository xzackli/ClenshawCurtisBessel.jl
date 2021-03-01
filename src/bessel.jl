

"""
    standard_bessel_moment(λ, ν, a)

Computes the standard Bessel moment ``G_{\\lambda, \\nu}(a)`` defined as
```math
G_{\\lambda, \\nu}(a) = \\int_0^a t^{\\lambda} J_{\\nu}(t)\\,dt.
```
"""
standard_bessel_moment(λ, ν, a) = pFq(SA[(1 + λ + ν)/2],SA[(3+λ+ν)/2, 1+ν], -a^2/4) * exp(
		log(a) * (1 + λ + ν) - log(2) * (1 + ν) +
		loggamma((1+λ+ν)/2) - loggamma((3+λ+ν)/2) - loggamma(1+ν))
