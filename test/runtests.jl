using AnalyticDiffEq
using Test

@testset "ExactDiffEq.jl" begin
    L(::Val{N}) where N<:Int = Ls[N-1]
    L0(x) = 0
    L1(x) = -x + 1
    L2(x) = .5*(x^2 -4x + 2)

    Ls = [L0 L1 L2]

    y = Variable()
    x = Variable()
    n = Variable()

    termlaguerre = @term(x*y'' + (1-x)*y' + n*y = 0)
    t1 = @term(x*y'')
    t2 = @term((1-x)*y')
    t3(i) = @term(n*y)


end
