module AnalyticDiffEq

using ModelingToolkit

abstract type AbstractEquation end
abstract type LaguerreAbstractEquation <: AbstractEquation end
abstract type HermiteAbstractEquation <: AbstractEquation end
abstract type LegendreAbstractEquation <: AbstractEquation end

struct LaguerreEquation{α,n} <: LaguerreAbstractEquation end
struct HermiteEquation{n} <: HermiteAbstractEquation end
struct LegendreEquation{n} <: LegendreAbstractEquation end

# struct HermiteEquation{n} <: AbstractEquation end
# struct LegendreEquation{n} <: AbstractEquation end
# struct ChebyshevFirstKindEquation{n} <: AbstractEquation end
# struct ChebyshevSecondKindEquation{n} <: AbstractEquation end
# struct GegenbauerEquation{α,n} <: AbstractEquation end
# struct JacobiEquation{α, β, n} <: AbstractEquation end

# Laguerre
# xy'' + (α + 1 - x)* y' + ny = 0
# eqtest1 = x*D2(y) + (1 + 1 - x) * D(y) + 1*y ~ 0
function isit(eq :: Equation, y :: Type{T}) where T<:AbstractEquation
        check_lhs(eq, y) && check_rhs(eq, y)
end

check_rhs(eq :: Equation, y :: Type{T}) where T<:AbstractEquation = eq.rhs.value == 0
check_lhs(eq :: Equation, y :: Type{T}) where T<:AbstractEquation = check_term1(eq,y) && check_term2(eq,y) && check_term3(eq,y)

function check_term1(eq :: Equation, y :: Type{LaguerreAbstractEquation})
        term1 = eq.lhs.args[1].args[1]
        flags = Bool[]

        flag1 = isequal(term1.op, *)
        push!(flags, flag1)

        flag2 = length(term1.args) == 2
        push!(flags, flag2)

        flag3 = typeof(term1.args[1]) == Variable
        push!(flags, flag3)

        x = term1.args[1]
        flag4 = x.known == false
        push!(flags, flag4)

        flag5 = x.name == :x # cross check dependents
        push!(flags, flag5)

        flag6 = length(x.dependents) == 1
        push!(flags, flag6)

        DD2 = term1.args[2]
        flag7 = typeof(DD2.op) == Differential
        push!(flags, flag7)

        # Time as a dependent?
        y = DD2.args[1].args[1]
        flag8 = isequal(y.dependents, x.dependents)
        push!(flags, flag8)

        flag9 = typeof(DD2.args[1].op) == Differential
        push!(flags, flag9)

        flag10 = y.name == :y
        push!(flags, flag10)

        flag11 = length(y.dependents) >= 1
        push!(flags, flag11)

        return all(flags)
end
function check_term2(eq :: Equation, y :: Type{LaguerreAbstractEquation})
        term2 = eq.lhs.args[1].args[2]
        flags = Bool[]

        flag1 = length(term2.args) == 2
        push!(flags, flag1)

        # check that it is indeed the derivative we want
        flag2 = typeof(term2.args[2].op) == Differential
        push!(flags, flag2)

        y = term2.args[2].args[1] # check names and dependents
        flag3 = length(y.dependents) == 1

        flag4 = typeof(y) == Variable
        push!(flags,flag4)

        # check it is an algebraic form we want
        flag5 = length(term2.args[1].args) == 2
        push!(flags, flag5)

        x = term2.args[1].args[2] # check names and dependents
        flag6 = typeof(x) == Variable
        push!(flags, flag6)

        # alpha + 1 TODO assume alpha is more than 0
        a = term2.args[1].args[1].value
        α = a - 1
        flag7 = α isa Real
        push!(flags, flag7)

        flag8 = term2.args[1].op == -
        push!(flags, flag8)

        return all(flags)
end
function check_term3(eq :: Equation, y :: Type{LaguerreAbstractEquation})
        term3 = eq.lhs.args[2]
        flags = Bool[]

        flag1 = isequal(term3.op, *)
        push!(flags, flag1)

        flag2 = length(term3.args) == 2
        push!(flags, flag2)

        y = term3.args[2]
        # TODO add check across terms so that all dependents are equal, all names are equal
        flag3 = typeof(y) == Variable
        push!(flags, flag3)

        flag4 = length(y.dependents) >= 1
        push!(flags, flag4)

        # Add logic for order 0
        n = term3.args[1].value
        flag5 = n isa Int && n >= 1
        push!(flags, flag5)

        return all(flags)
end

# Todo - cross check dependent variable is always the same
# This lets us use different dependent variabls as long as they are consistent
function check_dependents(eq :: Equation, y ::Type{LaguerreAbstractEquation})

end

# todo - cross check dependent names is always the same
function check_names(eq :: Equation, y :: Type{LaguerreAbstractEquation})

end


# Hermite
# y'' - 2xy' + λy = 0
# eq = D2(y) - 2*x*D(y) + λ*y ~ 0
function check_term1(eq :: Equation, y :: Type{HermiteAbstractEquation})
        term1 = eq.lhs.args[1].args[1]
        flags = Bool[]

        DD2 = term1
        flag1 = typeof(DD2.op) == Differential
        push!(flags, flag1)

        flag2 = typeof(DD2.args[1].op) == Differential
        push!(flags, flag2)

        y = term1.args[1].args[1]
        flag3 = length(y.dependents) >= 1
        push!(flags, flag3)

        flag4 = typeof(y) == Variable
        push!(flags, flag4)

        return all(flags)
end
function check_term2(eq :: Equation, y :: Type{HermiteAbstractEquation})
        term2 = eq.lhs.args[1].args[2]
        flags = Bool[]

        flag1 = length(term2.args) == 2
        push!(flags, flag1)

        # check that it is indeed the derivative we want
        flag2 = typeof(term2.args[2].op) == Differential
        push!(flags, flag2)

        y = term2.args[2].args[1] # check names and dependents
        flag3 = length(y.dependents) == 1

        flag4 = typeof(y) == Variable
        push!(flags,flag4)

        # check it is an algebraic form we want
        flag5 = length(term2.args[1].args) == 2
        push!(flags, flag5)

        x = term2.args[1].args[2] # check names and dependents
        flag6 = typeof(x) == Variable
        push!(flags, flag6)

        # alpha + 1 TODO assume alpha is more than 0
        flag7 = term2.args[1].args[1].value == 2
        push!(flags, flag7)

        flag8 = eq.lhs.args[1].op == -
        push!(flags, flag8)

        return all(flags)
end
function check_term3(eq :: Equation, y :: Type{HermiteAbstractEquation})
        term3 = eq.lhs.args[2]
        flags = Bool[]

        flag1 = isequal(term3.op, *)
        push!(flags, flag1)

        flag2 = length(term3.args) == 2
        push!(flags, flag2)

        y = term3.args[2]
        # TODO add check across terms so that all dependents are equal, all names are equal
        flag3 = typeof(y) == Variable
        push!(flags, flag3)

        flag4 = length(y.dependents) >= 1
        push!(flags, flag4)

        # Add logic for order 0
        n = term3.args[1].value
        flag5 = n isa Int && n >= 1
        push!(flags, flag5)

        return all(flags)
end

# TODO - cross check dependent variable is always the same
# This lets us use different dependent variabls as long as they are consistent
function check_dependents(eq :: Equation, y ::Type{HermiteAbstractEquation})

end

# TODO - cross check dependent names is always the same
function check_names(eq :: Equation, y :: Type{HermiteAbstractEquation})

end


# Legendre
# (1 - x^2) * y'' -2xy' + l(l+1)y = 0
# eqtest1 = (1-x^2) * D2(y) - 2*x*D(y) + l*(l+1)*y ~ 0
function check_term1(eq :: Equation, y :: Type{LegendreAbstractEquation})
        term1 = eq.lhs.args[1].args[1]
        flags = Bool[]

        flag1 = isequal(term1.op, *)
        push!(flags, flag1)

        flag2 = length(term1.args) == 2
        push!(flags, flag2)

        flag3 = typeof(term1.args[1]) == Operation
        push!(flags, flag3)

        parens = term1.args[1]
        flag4 = isequal(parens.op, -)
        push!(flags, flag4)

        flag5 = parens.args[1].value == 1
        push!(flags, flag5)

        xsq = parens.args[2]
        flag6 = length(xsq.args) == 2
        push!(flags, flag6)

        flag7 = isequal(xsq.op, ^)
        push!(flags, flag7)

        flag8 = xsq.args[2].value == 2
        push!(flags, flag8)

        DD2 = term1.args[2]
        flag9 = typeof(DD2.op) == Differential
        push!(flags, flag9)

        # Time as a dependent?
        y = DD2.args[1].args[1]
        flag10 = typeof(DD2.args[1].op) == Differential
        push!(flags, flag10)

        flag11 = length(y.dependents) >= 1
        push!(flags, flag11)

        return all(flags)
end
function check_term2(eq :: Equation, y :: Type{LegendreAbstractEquation})
        term2 = eq.lhs.args[1].args[2]
        flags = Bool[]

        flag1 = length(term2.args) == 2
        push!(flags, flag1)

        # check that it is indeed the derivative we want
        flag2 = typeof(term2.args[2].op) == Differential
        push!(flags, flag2)

        y = term2.args[2].args[1] # check names and dependents
        flag3 = length(y.dependents) == 1

        flag4 = typeof(y) == Variable
        push!(flags,flag4)

        # check it is an algebraic form we want
        flag5 = length(term2.args[1].args) == 2
        push!(flags, flag5)

        x = term2.args[1].args[2] # check names and dependents
        flag6 = typeof(x) == Variable
        push!(flags, flag6)

        flag8 = eq.lhs.args[1].op == -
        push!(flags, flag8)

        return all(flags)
end
function check_term3(eq :: Equation, y :: Type{LegendreAbstractEquation})
        term3 = eq.lhs.args[2]
        flags = Bool[]

        flag1 = isequal(term3.op, *)
        push!(flags, flag1)

        flag2 = length(term3.args) == 2
        push!(flags, flag2)

        y = term3.args[2]
        # TODO add check across terms so that all dependents are equal, all names are equal
        flag3 = typeof(y) == Variable
        push!(flags, flag3)

        flag4 = length(y.dependents) >= 1
        push!(flags, flag4)

        # Add logic for order 0
        l = term3.args[1].value
        flag5 = l isa Int && l >= 1 && llplus1checker(l)
        push!(flags, flag5)

        return all(flags)
end

function llplus1checker(l)
        l1 = isqrt(l)
        l2 = l1 + 1
        l1 * l2 == l ? true : false
end

# TODO - cross check dependent variable is always the same
# This lets us use different dependent variabls as long as they are consistent
function check_dependents(eq :: Equation, y ::Type{LegendreAbstractEquation})

end

# TODO - cross check dependent names is always the same
function check_names(eq :: Equation, y :: Type{LegendreAbstractEquation})

end


# Analytic Extraction

# Laguerre
# xy'' + (α + 1 - x)* y' + ny = 0
# eqlag = x*D2(y) + (1 + 1 - x) * D(y) + 1*y ~ 0
function extract(eq :: Equation, y :: Type{LaguerreAbstractEquation})
        n = eq.lhs.args[2].args[1].value
        α = eq.lhs.args[1].args[2].args[1].args[1].value - 1
        LaguerreEquation{α, n}()
end

# Hermite
# y'' - 2xy' + λy = 0
# eqher = D2(y) - 2*x*D(y) + λ*y ~ 0
function extract(eq :: Equation, y :: Type{HermiteAbstractEquation})
        n = eq.lhs.args[2].args[1].value
        HermiteEquation{n}()
end

# Legendre
# (1 - x^2) * y'' -2xy' + l(l+1)y = 0
# eqleg = (1-x^2) * D2(y) - 2*x*D(y) + l*(l+1)*y ~ 0
function extract(eq :: Equation, y :: Type{LegendreAbstractEquation})
        lp = eq.lhs.args[2].args[1].value
        n = isqrt(lp)
        LegendreEquation{n}()
end

function build(x :: LaguerreEquation{α,n}) where {α,n}
        return Ls[n + 1]
end

function build(x :: HermiteEquation{n}) where n
        return Hs[n + 1]
end

function build(x :: LegendreEquation{n}) where n
        Ps[n + 1]
end

function transform(x :: Equation, y ::Type{T}) where T<:AbstractEquation
        typ = extract(x,y)
        build(typ)
end

function analyze(eq :: Equation)
        try
                if isit(eq, LaguerreAbstractEquation)
                        return transform(eq, LaguerreAbstractEquation)
                end
        catch
        end

        try
                if isit(eq, HermiteAbstractEquation)
                        return transform(eq, HermiteAbstractEquation)
                end
        catch
        end

        try
                if isit(eq, LegendreAbstractEquation)
                        return transform(eq, LegendreAbstractEquation)
                end
        catch
        end

        throw("No analytic solution detected.")
end



export isit, extract, build, transform, analyze
export LaguerreAbstractEquation, LaguerreEquation
export HermiteAbstractEquation, HermiteEquation
export LegendreAbstractEquation, LegendreEquation

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

export Ls, Hs, Ps
end # module
