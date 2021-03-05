

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
M₃(a, ν) = ((32/a^3) * bessel_moment_G(3, ν, a) -
    (48/a^2) * bessel_moment_G(2, ν, a) +
    (18/a) * bessel_moment_G(1, ν, a) - bessel_moment_G(0, ν, a)) / a


@doc raw"""
    Mk_asymptotic(a::T, ν, k) where T

Computes a modified moment using asymptotics as ``k \rightarrow \infty``.
"""
function Mk_asymptotic(a::T, ν, k) where T

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
    BrandersPiessensProblem(a, ν)

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

Note that the coefficient indices decrease ``P_8, P_7, \cdots, P_0`` as the ``M_k``
index increases.
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
OliverR(OP::BrandersPiessensProblem{T}, i) where T = zero(T)


"""
    generate_BC(OP::BrandersPiessensProblem{T,8,2}, kmax::Int) where T

Generates the boundary conditions required for Oliver's method to solve
the boundary value problem posed Branders & Piessens.

# Arguments:
- `OP::BrandersPiessensProblem{T,8,2}`: problem parameters
- `kmax`: maximum order ``k`` to compute.

# Returns:
- `Tuple{Int,OffsetVector}`: the maximum index, and a sparse offset array
     with the boundaries.
"""
function generate_BC(OP::BrandersPiessensProblem{T,8,2}, kmax::Int) where T
    a, ν = OP.a, OP.ν
    ℓ = max(kmax+1, Int(ceil(a + 10)))  # added +1, never use asymptotic for kmax
    Mbc = OffsetVector(spzeros(T,ℓ+5), -3:(ℓ+1))

    Mbc[-3] = M₃(a, ν)
    Mbc[-2] = M₂(a, ν)
    Mbc[-1] = M₁(a, ν)
    Mbc[ 0] = M₀(a, ν)
    Mbc[ 1] = M₁(a, ν)
    Mbc[ 2] = M₂(a, ν)
    Mbc[ 3] = M₃(a, ν)

    Mbc[ℓ] = Mk_asymptotic(a, ν, ℓ)
    Mbc[ℓ+1] = Mk_asymptotic(a, ν, ℓ+1)

    return ℓ+1, Mbc
end



function momentM(Tout::Type{<:Number}, a::T, nu::T, kmax::Int, kbuffer=100) where T
    bpp = BrandersPiessensProblem{T,8,2}(a, nu)
    upperindex, Mbc = generate_BC(bpp, kmax+kbuffer)
    ρ = assembleρ(bpp, -3, upperindex, Mbc)
    𝐏 = assembleP(bpp, -3, upperindex)

    sol = OffsetVector(Array{Tout}(undef,kmax+1), 0:kmax)
    for i in 0:3
        sol[i] = Mbc[i]
    end
    res = qr(𝐏) \ ρ
    for i in 4:kmax
        sol[i] = Tout(res[i-3])
    end
    return sol
end
