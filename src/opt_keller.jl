
# currently unused, in principle could be used for computing higher order 
# Keller's method for integrating Bessel integrals through Chebyshev transform

@inline @fastmath A0(k, α::T, c::T, ν::T) where T = 2k * ( 2k^2 - 1 - 2α * (α+1) - c^2 + 2ν^2 )
@inline @fastmath A1(k, α::T, c::T, ν::T) where T = 2k^3 - (1 + 4α) * k^2 + (2α + 2α^2 - c^2 / 2 - 2ν^2) * k - α^2 - 5c^2/4
@inline @fastmath A2(k, α::T, c::T, ν::T) where T = c^2 * (k-1)
@inline @fastmath A3(k, α::T, c::T, ν::T) where T = T(1/4) * c^2 * (2k - 1)


@inline Ak_coeffs(k, α::T, c::T, ν::T) where T = SVector{7,T}(
    -A3(-k, α, c, ν), -A2(-k, α, c, ν), -A1(-k, α, c, ν), A0(k, α, c, ν), 
    A1(k, α, c, ν), A2(k, α, c, ν), A3(k, α, c, ν))

Base.@propagate_inbounds @fastmath r_rhs(k::Int, kf, bₖ) = 4kf * bₖ[k] - (2kf+1) * bₖ[k-1] - (2kf-1) * bₖ[k+1]


Base.@propagate_inbounds function _ge_kernel_CASE_A(f::T, L1::SVector{7,T},  L2::SVector{7,T}) where T
    return @fastmath SVector{7,T}(
        muladd(f, L1[2], L2[2]), 
        muladd(f, L1[3], L2[3]),
        muladd(f, L1[4], L2[4]),
        muladd(f, L1[5], L2[5]), 
        muladd(f, L1[6], L2[6]), 
        muladd(f, L1[7], L2[7]), zero(T))
end

# case A: m = 0, s = -1
function ge_CASE_A(L::Vector{SVector{7,T}}, b::SVector, r::Vector{T}, f::Vector{T}, 
                   α::T, c::T, ν::T, ϵ::T; min_stop=3) where T
	m = length(L)
    @assert m > 3

    k₁, k₂, k₃ = T(3), T(4), T(5)
    ℵ₁ = SVector{7,T}(-A1(-k₁, α, c, ν), A0(k₁, α, c, ν), A1(k₁, α, c, ν), A2(k₁, α, c, ν), A3(k₁, α, c, ν), zero(T), zero(T))
    ℵ₂ = SVector{7,T}(-A2(-k₂, α, c, ν), -A1(-k₂, α, c, ν), A0(k₂, α, c, ν), A1(k₂, α, c, ν), A2(k₂, α, c, ν), A3(k₂, α, c, ν), zero(T))
    ℵ₃ =  Ak_coeffs(k₃, α, c, ν)
    ϵ̃ = ϵ * T(1/30)  # stopping condition, multiplying the 1/30 out of the loop

    @inbounds for i ∈ 1:(m-3)  # set up RHS
        r[i] = r_rhs(i+2, T(i+2), b)  # first is k=3
    end

    @inbounds for i ∈ 1:(m-3)

        # check if we need to pivot
        k₁, k₂, k₃ = T(i+2), T(i+3), T(i+4)  # start from k = 3
        first_ℵ₁, first_ℵ₂, first_ℵ₃ = first(ℵ₁), first(ℵ₂), first(ℵ₃)
        abs_first_ℵ₁, abs_first_ℵ₂, abs_first_ℵ₃ = abs(first_ℵ₁), abs(first_ℵ₂), abs(first_ℵ₃)

        if abs_first_ℵ₂ > abs_first_ℵ₁ && abs_first_ℵ₂ > abs_first_ℵ₃        # have to pivot on second
            ℵ₁, ℵ₂ = ℵ₂, ℵ₁
            first_ℵ₁, first_ℵ₂ = first_ℵ₂, first_ℵ₁
            abs_first_ℵ₁, abs_first_ℵ₂ = abs_first_ℵ₂, abs_first_ℵ₁
            r[i], r[i+1] = r[i+1], r[i]
        elseif abs_first_ℵ₃ > abs_first_ℵ₁ && abs_first_ℵ₃ > abs_first_ℵ₂    # have to pivot on third
            ℵ₁, ℵ₃ = ℵ₃, ℵ₁
            first_ℵ₁, first_ℵ₃ = first_ℵ₃, first_ℵ₁
            abs_first_ℵ₁, abs_first_ℵ₃ = abs_first_ℵ₃, abs_first_ℵ₁
            r[i], r[i+2] = r[i+2], r[i]
        end
        
        # gaussian elimination after pivoting
        L[i] = ℵ₁
        p = T(-1) / first_ℵ₁  # negative inverse of pivot element, elimination is ADDITION
        f[i] = -p  # store negative reciprocals in f for the backsolve later
        f₂ = p * first_ℵ₂
        f₃ = p * first_ℵ₃
        ℵ₂ = _ge_kernel_CASE_A(f₂, ℵ₁, ℵ₂)
        ℵ₃ = _ge_kernel_CASE_A(f₃, ℵ₁, ℵ₃)

        # also row reduce the rhs
        r₁ = r[i]
        r₂ = r[i+1] + f₂ * r₁
        r₃ = r[i+2] + f₃ * r₁
        r[i+1] = r₂
        r[i+2] = r₃

        # check if we should stop
        if i > min_stop
            r₁, r₂, r₃ = abs(r₁), abs(r₂), abs(r₃)
            end_cond = k₁^2 * (r₁ + r₂ + r₃) # max(r₁, r₂, r₃)
            if end_cond < ϵ̃ * abs_first_ℵ₁
                return i
            end
        end

        # set up for next time
        ℵ₁, ℵ₂ = ℵ₂, ℵ₃
        ℵ₃ = Ak_coeffs(T(i+5), α, c, ν)
    end

    return (m-2)
end

# we store the reciprocal diagonal in f when performing Gaussian elimination
Base.@propagate_inbounds @fastmath function ge_backsolve!(n, L, r, f)
    f[n] *= r[n]
    for i in (n-1):-1:1
        yᵢ = r[i]
        Lᵢ = L[i]
        for Δ ∈ 1:min(6, n - i)  # banded, at most 7 elts
            yᵢ -= Lᵢ[1 + Δ] * f[i + Δ]
        end
        f[i] *= yᵢ
    end
    f
end


# Base.@propagate_inbounds @fastmath function ge_backsolve!_old(n, L, r, f)
#     f[n] = r[n] / L[n][1]
#     for i in (n-1):-1:1
#         yᵢ = r[i]
#         Lᵢ = L[i]
#         for Δ ∈ 1:min(6, n - i)  # banded, at most 7 elts
#             yᵢ -= Lᵢ[1 + Δ] * f[i + Δ]
#         end
#         f[i] = yᵢ / Lᵢ[1]
#     end
#     f
# end


# debug test routines
function debug_Aswitch(s, k, α, c, ν)
    if s == -3
        return -A3(-k, α, c, ν)
    elseif s == -2
        return -A2(-k, α, c, ν)
    elseif s == -1
        return -A1(-k, α, c, ν)
    elseif s == 0
        return A0(k, α, c, ν)
    elseif s == 1
        return A1(k, α, c, ν)
    elseif s == 2
        return A2(k, α, c, ν)
    elseif s == 3
        return A3(k, α, c, ν)
    else
        return 0.
    end
end

function debug_materialize_A(n, α, c, ν)
    A = zeros(n, n)
    for i in 1:n, j in 1:n
        A[i, j] = debug_Aswitch(j - i - 1, i + 2, α, c, ν)
    end
    A
end

function debug_materialize_r(n, b)
    r = zeros(n)
    for i in 1:(n-4)
        r[i] = r_rhs(i+2, Float64(i+2), b)  # first is k=3
    end
    r
end

