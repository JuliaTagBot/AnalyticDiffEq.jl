using Test
using ModelingToolkit
using AnalyticDiffEq
# Laguerre
# xy'' + (α + 1 - x)* y' + ny = 0

# Hermite
# y'' - 2xy' + λy = 0

# Legendre
# (1 - x^2) * y'' -2xy' + l(l+1)y = 0

L0(x) = 1.0
    L1(x) = -x + 1.0
    L2(x) = .5 * (x^2 - 4x + 2)
    L3(x) = (1/6)*(-x^3 +9x^2 -18x + 6)
    L4(x) = (1/24)*(x^4 - 16x^3 + 72x^2 - 96x + 24)
    L5(x) = (1/120)*(-x^5 + 25x^4 -200x^3 + 600x^2 - 600x + 120)
    L6(x) = (1/720)*(x^6 - 36x^5 + 450x^4 - 2400x^3 + 5400x^2 - 4320x + 720)
    L7(x) = (1/5040)*(-x^7 + 49x^6 - 882x^5 + 7350x^4 - 29_400x^3 + 52_920x^2 - 35_280x + 5040)
    L8(x) = (1/40320)*(x^8 - 64x^7 + 1568x^6 - 18_816x^5 + 117_600x^4 - 376_320x^3 + 564_480x^2 - 322_560x + 40320)
    L9(x) = (1/362880)*(-x^9 + 81x^8 - 2592x^7 + 42_336x^6 - 381_024x^5 + 1_905_120x^4 - 5_080_320x^3 + 6_531_840x^2 - 3_265_920x + 362_880)
    L10(x) = (1/3628800)*(x^10 - 100x^9 + 4_050x^8 - 86_400x^7 + 1_058_400x^6 - 7_620_480x^5 + 31_752_000x^4 - 72_576_000x^3  + 81_648_000x^2 - 36_288_000x + 3_628_800)
    L11(x) = (1/39916800)*(-x^11 + 121x^10 - 6_050x^9 + 163_350x^8 - 2_613_600x^7 + 25_613_280x^6 - 153_679_680x^5 + 548_856_000x^4 - 1_097_712_000x^3 + 1_097_712_000x^2 - 439_084_800x + 39_916_800)
    L12(x) = (1/479001600)*(x^12 - 144x^11 + 8_712x^10 - 290_400x^9 + 5_880_600x^8 - 75_271_680x^7 + 614_718_720x^6 - 3_161_410_560x^5 + 9_879_408_000x^4 - 17_563_392_000x^3 + 15_807_052_800x^2 - 5_748_019_200x + 479_001_600)

    Ls = [L0, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12]

    # Credit to @yingboma
    function laguerre(n, x)

        n == 0 ? (return 1) : nothing

        l0, l1 = 1, 1-x
        for k in 1:n-1
            l0, l1 = l1, ( (2k+1-x)*l1 - k*l0 )/(k+1)
        end
        l1
    end

H0(x) = 1.0
    H1(x) = 2x
    H2(x) = -2 + 4x^2
    H3(x) = -12x + 8x^3
    H4(x) = 12 - 48x^2 + 16x^4
    H5(x) = 120x -160x^3 + 32x^5
    H6(x) = -120 + 720x^2 - 480x^4 + 64x^6
    H7(x) = -1680x + 3360x^3 - 1344x^5 + 128x^7
    H8(x) =1680 - 13440x^2 + 13440x^4 - 3584x^6 + 256x^8
    H9(x) =30240x - 80640x^3 + 48384x^5 - 9216x^7 + 512x^9
    H10(x) =-30240 + 302400x^2 - 403200x^4 + 161280x^6 - 23040x^8 + 1024x^10
    H11(x) =-665280x + 2217600x^3 - 1774080x^5 + 506880x^7 - 56320x^9 + 2048x^11
    H12(x) =665280 - 7983360x^2 + 13305600x^4 - 7096320x^6 + 1520640x^8 - 135168x^10 + 4096x^12

    Hs = [H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12]

    function hermite(x,n)
        if n == 0
            return 1
        elseif n == 1
            return 2x
        else
            return 2x*hermite(x,n-1) - 2*(n-1)*hermite(x,n-2)
        end
    end

P0(x) = 1
    P1(x) =  x
    P2(x) = .5*(3x^2 -1)
    P3(x) = .5*(5x^3 -3x)
    P4(x) = (1/8)*(35x^4 - 30x^2 + 3)
    P5(x) = (1/8)*(63x^5 - 70x^3 + 15x)
    P6(x) = (1/16)*(231x^6 - 315x^4 + 105x^2 - 5)
    P7(x) = (1/16)*(429x^7 - 693x^5 + 315x^3 - 35x)
    P8(x) = (1/128)*(6435x^8 - 12012x^6 + 6930x^4 - 1260x^2 + 35)
    P9(x) = (1/128)*(12155x^9 - 25740x^7 + 18018x^5 - 4620x^3 + 315x)
    P10(x) = (1/256)*(46189x^10 - 109395x^8 + 90090x^6- 30030x^4 + 3465x^2 - 63)
    P11(x) = (1/256)*(2^11)*(705432x^11 - 1847560x^9 + 1750320x^7 - 720720x^5 + 120120x^3 - 5544x)
    P12(x) = (1/512)*(2^12)*(2704156x^12 - 7759752x^10 + 8314020x^8 - 4084080x^6 + 900900x^4 - 72072x^2 + 924)

    Ps = [P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12]
    testPs = [P2, P3, P9, P10]

    function legendre(x, n)
        if n == 0
            return 1
        elseif n == 1
            return x
        else
            return ((2n-1)*x*legendre(x,n-1) - (n - 1)*legendre(x, n-2))/n
        end
    end

@testset "AnalyticDiffEq.jl" begin

    @testset "Analytic Detection" begin

        # Laguerre
        # xy'' + (α + 1 - x)* y' + ny = 0
        # eq = x*D2(y) + (α + 1 - x) * D(y) + n*y ~ 0
        @testset "Laguerre" begin

            @parameters t
            @variables y(t) x(t)
            @derivatives D'~t D2''~t

            eqtest1 = x*D2(y) + (1 + 1 - x) * D(y) + 1*y ~ 0
            eqtest2 = x*D2(y) + (6 - x) * D(y) + 3*y ~ 0
            eqtest3 = x*D2(y) + (1 - x) * D(y) + 7*y ~ 0
            # Weird corner cases
            # eqtest4 = x*D2(y) + (7 - x) * D(y) ~ 0
            # eqtest0 = x*D2(y) + (1 - x) * D(y) ~ 0
            # eqtest5 = x*D2(y) + (1 - x) * D(y) + y ~ 0

            LAE = LaguerreAbstractEquation

            @test isit(eqtest1, LAE)
            @test isit(eqtest2, LAE)
            @test isit(eqtest3, LAE)
            # @test_broken isit(eqtest4, LaguerreAbstractEquation)
            # @test isit(eqtest5, LaguerreAbstractEquation)
            # @test_broken isit(eqtest0, LaguerreAbstractEquation)

       end

       # Hermite
       # y'' - 2xy' + λy = 0
       # eqtest1 = D2(y) - 2*x*D(y) + λ*y ~ 0
        @testset "Hermite" begin

            @parameters t
            @derivatives D'~t D2''~t
            @variables y(t) x(t)

            eqtest1 = D2(y) - 2*x*D(y) + 1*y ~ 0
            eqtest2 = D2(y) - 2*x*D(y) + 7*y ~ 0
            eqtest3 = D2(y) - 2*x*D(y) + 12*y ~ 0

            HAE = HermiteAbstractEquation

            @test isit(eqtest1, HAE)
            @test isit(eqtest2, HAE)
            @test isit(eqtest3, HAE)

        end

        # Legendre
        # (1 - x^2) * y'' -2xy' + l(l+1)y = 0
        # eqtest1 = (1-x^2) * D2(y) - 2*x*D(y) + l*(l+1)*y ~ 0
        @testset "Legendre" begin

            @parameters t
            @derivatives D'~t D2''~t
            @variables y(t) x(t)

            eqtest1 = (1-x^2) * D2(y) - 2*x*D(y) + 2*y ~ 0
            eqtest2 = (1-x^2) * D2(y) - 2*x*D(y) + 6*y ~ 0
            eqtest3 = (1-x^2) * D2(y) - 2*x*D(y) + 12*y ~ 0

            @test AnalyticDiffEq.llplus1checker(6) == true
            @test AnalyticDiffEq.llplus1checker(2) == true
            @test AnalyticDiffEq.llplus1checker(12) == true
            @test AnalyticDiffEq.llplus1checker(20) == true

            LegAE = LegendreAbstractEquation
            @test isit(eqtest1, LegAE)
            @test isit(eqtest2, LegAE)
            @test isit(eqtest3, LegAE)

        end
    end

    @testset "Analytic Extraction" begin

        @parameters t
        @derivatives D'~t D2''~t
        @variables y(t) x(t)

        eqlag = x*D2(y) + (0 + 1 - x) * D(y) + 6*y ~ 0
        eqher = D2(y) - 2*x*D(y) + 5*y ~ 0
        eqleg = (1-x^2) * D2(y) - 2*x*D(y) + 4*(4+1)*y ~ 0
        @test extract(eqlag, LaguerreAbstractEquation) == LaguerreEquation{0,6}()
        @test extract(eqher, HermiteAbstractEquation) == HermiteEquation{5}()
        @test extract(eqleg, LegendreAbstractEquation) == LegendreEquation{4}()

    end

    @testset "Analytic Construction" begin

        @parameters t
        @derivatives D'~t D2''~t
        @variables y(t) x(t)

        eqlag = x*D2(y) + (0 + 1 - x) * D(y) + 6*y ~ 0
        eqher = D2(y) - 2*x*D(y) + 5*y ~ 0
        eqleg = (1-x^2) * D2(y) - 2*x*D(y) + 4*(4+1)*y ~ 0
        eqerr = 0 ~ x^2
        LagEQ = extract(eqlag, LaguerreAbstractEquation)
        HerEQ = extract(eqher, HermiteAbstractEquation)
        LegEQ = extract(eqleg, LegendreAbstractEquation)

        @test build(LagEQ)(.3) ≈ Ls[7](.3)
        @test build(HerEQ)(.3) ≈ Hs[6](.3)
        @test build(LegEQ)(.3) ≈ Ps[5](.3)

        @test analyze(eqlag)(.3) ≈ Ls[7](.3)
        @test analyze(eqher)(.3) ≈ Hs[6](.3)
        @test analyze(eqleg)(.3) ≈ Ps[5](.3)
        @test_broken analyze(eqerr) == 3.14

    end

    @testset "Solve Interface" begin
        @test_broken 1 == 2
    end
end
