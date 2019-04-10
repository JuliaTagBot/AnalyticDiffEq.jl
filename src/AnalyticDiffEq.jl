module AnalyticDiffEq

using Terms

y = Variable()
x = Variable()
n = Variable()

termlaguerre = @term(x*y'' + (1-x)*y' + n*y = 0)
target = @term(x*y'' + (1-x)*y' + 3*y = 0)
tl.head

nlaguerre(t) = t.args[1].args[4].args[2].head

function ist1laguerre(t)
    t.args[1].args[2].args[2].head == Variable()

t1 = @term(x*y'')
t2 = @term((1-x)*y')
t3 = @term(n*y)

tl

tsl(t1,t2) = t1 ∈ children(t2.args[1]) ? true : false

ist1laguerre(t) = tsl(t1,t)
ist2laguerre(t) = tsl(t2,t)
ist3laguerre(t) =
ist4laguerre(t) = t.args[2].args[2].head == 0 ? true : false
isheadlaguerre(t) = t == :(=) ? true : false
termslaguerre(t) = ist1laguerre(t) && ist2laguerre(t) && ist3laguerre(t)
remainderlaguerre(t) = ist4laguerre(t) && isheadlaguerre(t)
isitlaguerre(t) = termslaguerre(t) && remainderlaguerre(t)

end # module
