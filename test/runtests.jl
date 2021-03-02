using ClenshawCurtisBessel
using Test
using DoubleFloats
using LinearAlgebra

# """
# Consider an example recurrence relation, which we will call `Example61161`
# ```math
# 6 y(s) - 11 y(s+1) + 6 y(s+2) - y(s+3) = 0.
# ```
# with initial conditions ``y(1) = 2``, ``y(2)=5``, and ``y(3) = 15``.
# """
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

    @test ClenshawCurtisBessel.Mk_asymptotic(2.0, 3.0, 21.0) ≈ -0.0001477774680739161 rtol=1e-6
    @test ClenshawCurtisBessel.Mk_asymptotic(3.0, 5.0, 21.0) ≈ -0.00004961332181983246 rtol=1e-6
    @test ClenshawCurtisBessel.Mk_asymptotic(2.0, 6.0, 24.0) ≈ -1.06156879586713e-6 rtol=1e-6
    @test ClenshawCurtisBessel.Mk_asymptotic(3., 7.0, 28.0) ≈ -1.647048323970303e-6 rtol=1e-6
end


@testset "Double64 Raw Moments" begin
	CCB = ClenshawCurtisBessel
    bpp = CCB.BrandersPiessensProblem{Double64,8,2}(10.0, 3 + 1/2)
	indexBC, MBC = CCB.generate_BC(bpp, 256)
	ρ = CCB.assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = CCB.assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ
	ref = big"0.0002069511037367724863632484164263304887654"
	@test (sol[11] .- ref) ./ ref < 1e-15

    bpp = CCB.BrandersPiessensProblem{Double64,8,2}(10.0, 3)
	indexBC, MBC = CCB.generate_BC(bpp, 256)
	ρ = CCB.assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = CCB.assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ
	ref = big"-0.0001955581634668262545942178152189286110836"
	@test (sol[11] .- ref) ./ ref < 1e-15
	ref = big"-0.00009457154862533141412967956613252113399756"
	@test (sol[16] .- ref) ./ ref < 1e-15
	ref = big"-0.00001631452917562497345315040327964024786504"
	@test (sol[40] .- ref) ./ ref < 1e-15
end

@testset "Float64 Raw Moments" begin
	CCB = ClenshawCurtisBessel
    bpp = CCB.BrandersPiessensProblem{Float64,8,2}(10.0, 3 + 1/2)
	indexBC, MBC = CCB.generate_BC(bpp, 256)
	ρ = CCB.assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = CCB.assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ
	ref = big"0.0002069511037367724863632484164263304887654"
	@test (sol[11] .- ref) ./ ref < 1e-6

    bpp = CCB.BrandersPiessensProblem{Float64,8,2}(10.0, 3)
	indexBC, MBC = CCB.generate_BC(bpp, 256)
	ρ = CCB.assembleρ(bpp, -3, indexBC, MBC)
	𝐏 = CCB.assembleP(bpp, -3, indexBC)
	sol = qr(𝐏) \ ρ
	ref = big"-0.0001955581634668262545942178152189286110836"
	@test (sol[11] .- ref) ./ ref < 1e-6
end


@testset "Double64 momentM" begin

    ref = [big"0.0007290114582988161923809468346226923322190",
        big"0.0005938052788264036365639278651356395999367",
        big"0.0002786292946361793775883850681999983130955",
        big"-0.00002805260432115378687369846500283467077757",
        big"-0.0001950344048967695886304586424720522836879",
        big"-0.0002181086699578411818278638336785129268495",
        big"-0.0001688344453897982834345306352230680868134",
        big"-0.0001124694826392090914811134620250485009294",
        big"-0.00007437704414078541354435912139888596042476",
        big"-0.00005298335002672124481290725819716655026908",
        big"-0.00004066037955955940010097890154609528857357",
        big"-0.00003259729530260456078582391639304253599224",
        big"-0.00002680486467601182798029350986861852376451",
        big"-0.00002246247645744408637896618989749724244207",
        big"-0.00001911995263563955270424343152423340974346",
        big"-0.00001648755827827632619531486527620890550541",
        big"-0.00001437360220319041810429287944102381603564",
        big"-0.00001264808824960577566084971381515466084373",
        big"-0.00001121999433319545122043778233304758581356",
        big"-0.00001002382210812218822542518716028642086529",
        big"-9.011366672874117063012554033327179984774e-6"]
    # note that res starts at 0
    res = momentM(Float64, DoubleFloat(6.0), DoubleFloat(10.0), 20)
    for i in 0:19
        @test abs(res[i] .- ref[i+1]) < 1e-16
    end

    # now test the half-integer case
    ref = [big"0.0003784769456217184219402128017246464332000",
        big"0.0003116266706253777562730807652926920232194",
        big"0.0001539041845785743007034880912831886326126",
        big"-3.947735731437409097894699973916750717989e-6",
        big"-0.00009548112619274286495691961673995490768193",
        big"-0.0001140573571560400900107740517746052585751",
        big"-0.00009226229577888985713288717817974834079267",
        big"-0.00006329298661640707880778541216716573215075",
        big"-0.00004220681212615568225179049644556672266289",
        big"-0.00002983588500153932889283750688859037769246",
        big"-0.00002268346477382750545447550291183506504071"]
    # note that res starts at 0
    res = momentM(Float64, DoubleFloat(6.0), DoubleFloat(10.5), 20)
    for i in 0:10
        @test abs(res[i] .- ref[i+1]) < 1e-16
    end
end
