using ClenshawCurtisBessel
using Test


@testset "Oliver's Method" begin
    y61161(s) = 1 - 2^(s-1) + 2 * 3^(s-1)
    struct Example61161{T,N,M} <: OliverProblem{T,N,M} end
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
    for a in 0:20
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
                P = assembleP(ex, a, b)
                ρ = assembleρ(ex, a, b, YBC)

                # """$(P * y.((a+n-m+1):(b-m))) $(ρ)"""
                sol = P \ ρ
                ref = y61161.((a+n-m+1):(b-m))
                @test all( ((sol .- ref) ./ ref) .< 1e-11 )
            end
        end
    end
end

@testset "Moments" begin
    @test ClenshawCurtisBessel.standard_bessel_moment(10, 15/2, 2) ≈ 0.007089503058490749
    @test ClenshawCurtisBessel.standard_bessel_moment(100, 15/2, 2) ≈ 1.482273470904133e24
    @test ClenshawCurtisBessel.standard_bessel_moment(70, 59/2, 1) ≈ 2.695679919047616e-43
    @test ClenshawCurtisBessel.standard_bessel_moment(100, 5000, 1) ≈ 0.0
    @test ClenshawCurtisBessel.standard_bessel_moment(1, 5, 1) ≈ 0.00003601420867202275
end
