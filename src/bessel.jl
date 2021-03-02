

"""
    bessel_moment_G(λ, ν, a)

Computes the standard Bessel moment ``G_{\\lambda, \\nu}(a)`` defined as
```math
G_{\\lambda, \\nu}(a) = \\int_0^a t^{\\lambda} J_{\\nu}(t)\\,dt.
```
"""
function bessel_moment_G(λ, ν, a::T) where T
    return HypergeometricFunctions.pFq(
            SA[(1 + λ + ν)/2],SA[(3+λ+ν)/2, 1+ν], (-a^2/4)) * exp(
        log(a) * (1 + λ + ν) - log(2) * (1 + ν) +
        lgamma((1+λ+ν)/2) - lgamma((3+λ+ν)/2) - lgamma((1+ν)))
end

# modified moment boundary conditions (Piessens & Branders 1983 eq. 15-18)
M₀(a, ν) = bessel_moment_G(0, ν, a) / a
M₁(a, ν) = ((2/a) * bessel_moment_G(1, ν, a) - bessel_moment_G(0, ν, a)) / a
M₂(a, ν) = ((8/a^2)*bessel_moment_G(2, ν, a) -
    (8/a)*bessel_moment_G(1, ν, a) + bessel_moment_G(0, ν, a)) / a
M₃(a, ν) = (
    (32/a^3) * bessel_moment_G(3, ν, a) -
    (48/a^2) * bessel_moment_G(2, ν, a) +
    (18/a) * bessel_moment_G(1, ν, a) - bessel_moment_G(0, ν, a)
) / a


@doc raw"""
    Mₖ_asymptotic(a::T, ν, k) where T

Computes a modified moment using asymptotics as ``k \rightarrow \infty``.
"""
function Mₖ_asymptotic(a::T, ν, k) where T

    ϕ0 = SA[
        J(ν, a),
        zero(T),
        ((3ν-2) * J(ν, a) - 3a * J(ν-1, a)) / 2,
        zero(T),
        ( 15a * J(ν-1, a) + (4 - 15a^2 + 15ν * (ν-1)) * J(ν, a) ) / 4,
        0,
        (21a * (-2 + 5a^2 - 5ν^2) * J(ν-1, a) +
            (-8 - 105a^2 * (ν-3) + 21ν * (2 + 5 * (ν-1) * ν)) * J(ν, a)) / 8
    ]

    ϕπ = SA[
        -J(ν, zero(T)),
        zero(T),
        J(ν, zero(T)) - (3a/4) * J(ν-1, zero(T)),
        zero(T),
        (-15a^2 * J(ν-2, zero(T)) + 60a * J(ν-1, zero(T)) +
            2 * (15a^2 - 8) * J(ν, zero(T))) / 16,
        zero(T),
        (-105a^3 * J(ν-3, zero(T)) +
            840a^2 * J(ν-2, zero(T)) -
            (315a^3 - 1008a) * J(ν-1, zero(T)) +
            (64 - 1680a^2) * J(ν, zero(T))
        ) / 64
    ]

    Σ = zero(T)
    for j in 0:3
        Σ = Σ + (-1)^j * ((-1)^k * ϕπ[2j+1] - ϕ0[2j+1]) * T(k)^(-2j-2)
    end

    return Σ / 2
end


"""
Construct a boundary value problem from the modified moment recurrence relation.

`N` and `M` correspond to the total number of specified initial boundary conditions,
and the number specified on the end of the recurrence, respectively. T is the type
we're computing on, `a` is the parameter in the Bessel function Jᵥ(ax), and `ν` is
the order.
"""
struct BrandersPiessensProblem{T,N,M} <: OliverProblem{T,N,M}
    a::T
    ν::T
end


@doc raw"""
    ClenshawCurtisBessel.OliverP(OP::BrandersPiessensProblem{T}, s, k) where T

Computes the modifed moments from Branders & Piessens 1987.

```math
\begin{aligned}
\frac{a^2}{16} &M_{k}(a, \nu) + \left[ (k+1)^2 - \nu^2 - \frac{a^2}{4} \right]
M_{k + 2}(a, \nu) + \left[ 4 \nu^2 - 2(k+4) + 4 \right] M_{k+3}(a,\nu) \\
&- \left[ 2(k+4)^2 - 6 + 6\nu^2 - \frac{3a^2}{8} \right] M_{k+4}(a,\nu) +
(4\nu^2 + 2(k+4) + 4)M_{k+5}(a,\nu) \\
&+ \left[ (k+7)^2 - \nu^2 - \frac{a^2}{4} \right]M_{k+6} + \frac{a^2}{16} M_{k+8}(a,\nu) \quad = 0
\end{aligned}
```
"""
function OliverP(
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

"""
    OliverR(OP::BrandersPiessensProblem{T}, i)

Computes ``R(i)``, the right side of the recurrence relation.
"""
OliverR(
    OP::BrandersPiessensProblem{T}, i) where T = zero(T)

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
    MBC[ℓ] = Mₖ_asymptotic(a, ν, ℓ)
    MBC[ℓ+1] = Mₖ_asymptotic(a, ν, ℓ+1)

    return ℓ+1, MBC
end
