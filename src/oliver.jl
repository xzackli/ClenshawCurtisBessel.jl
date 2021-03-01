
"""
Construct a boundary value problem from a recurrence relation. The functions for
coefficients and the right side dispatch on this abstract type.

`N` and `M` correspond to the total number of specified initial boundary conditions,
and the number specified on the end of the recurrence, respectively.
"""
abstract type OliverProblem{T,N,M} end

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

"""
    assembleρ{T}(OP::OliverProblem{M,N}, a, b) where {M, N}

Generic description of the ``\\vec{\rho}`` vector within the linear system of Oliver's method.
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


"""
	OliverP(OP::OliverProblem, s, i)

Computes ``P_s(i)``, the coefficient ``s`` at solution index ``i``.
"""
function OliverP(OP::OliverProblem, s, i)
    throw(ErrorException("You need to implement an OliverP for this problem."))
end

"""
	OliverR(OP::OliverProblem, i)

Computes ``R(i)``, the right side of the recurrence relation.
"""
function OliverR(OP::OliverProblem, i)
    throw(ErrorException("You need to implement an OliverR for this problem."))
end
