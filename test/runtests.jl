using ClenshawCurtisBessel
using Test

"""
Consider an example recurrence relation, which we will call `Example61161`
```math
6 y(s) - 11 y(s+1) + 6 y(s+2) - y(s+3) = 0.
```
with initial conditions ``y(1) = 2``, ``y(2)=5``, and ``y(3) = 15``.
"""
@testset "Oliver's Method" begin
    y61161(s) = 1 - 2^(s-1) + 2 * 3^(s-1)
    struct Example61161{T,N,M} <: ClenshawCurtisBessel.OliverProblem{T,N,M} end
    function ClenshawCurtisBessel.OliverP(OP::Example61161, s, i)
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
    ClenshawCurtisBessel.OliverR(OP::Example61161, i) = 0

    # loop over a broad range of BC and sizes
    for a in -2:10
        for b in (a+3):(a+400)
            for m in 0:3
                a, b, n = 1, 9, 3

                YBC = Dict{Int,Int}()
                for i in (a+1):(a+n-m)
                    YBC[i] = y61161(i)
                end
                for i in (b-m+1):b
                    YBC[i] = y61161(i)
                end

                ex = Example61161{Float64,n,m}()
                P = ClenshawCurtisBessel.assembleP(ex, a, b)
                ρ = ClenshawCurtisBessel.assembleρ(ex, a, b, YBC)

                # """$(P * y.((a+n-m+1):(b-m))) $(ρ)"""
                sol = P \ ρ
                ref = y61161.((a+n-m+1):(b-m))
                @test all( ((sol .- ref) ./ ref) .< 1e-11 )
            end
        end
    end
end

# MATHEMATICA INPUT FOR THESE INTEGRALS
# NumberForm[NIntegrate[1/2 ChebyshevT[k, x] BesselJ[\[Nu], (x + 1) a/2] /.
# {a -> 10, \[Nu] -> 3 , k -> 3 + 11} , {x, -1, 1}, AccuracyGoal -> 40,
# PrecisionGoal -> 40, WorkingPrecision -> 40], 40]
@testset "Moment Boundary Conditions" begin
    @test ClenshawCurtisBessel.bessel_moment_G(10, 15/2, 2) ≈ 0.007089503058490749
    @test ClenshawCurtisBessel.bessel_moment_G(100, 15/2, 2) ≈ 1.482273470904133e24
    @test ClenshawCurtisBessel.bessel_moment_G(70, 59/2, 1) ≈ 2.695679919047616e-43
    @test ClenshawCurtisBessel.bessel_moment_G(100, 5000, 1) ≈ 0.0
    @test ClenshawCurtisBessel.bessel_moment_G(1, 5, 1) ≈ 0.00003601420867202275

    @test ClenshawCurtisBessel.Mₖ_asymptotic(2.0, 3.0, 21.0) ≈ -0.0001477774680739161 rtol=1e-6
    @test ClenshawCurtisBessel.Mₖ_asymptotic(3.0, 5.0, 21.0) ≈ -0.00004961332181983246 rtol=1e-6
    @test ClenshawCurtisBessel.Mₖ_asymptotic(2.0, 6.0, 24.0) ≈ -1.06156879586713e-6 rtol=1e-6
    @test ClenshawCurtisBessel.Mₖ_asymptotic(3., 7.0, 28.0) ≈ -1.647048323970303e-6 rtol=1e-6
end


@testset "Moments" begin

end
