# AnalyticDiffEq - WIP!

### Why estimate when you can solve?

In some very special cases (usually covered in introductory differential equations textbooks), differential equations have a magical formula that can detect its analytic solution.

With the Julia programming language and the full power of the DifferentialEquations.jl ecosystem, we can build on top of ModelingToolkit.jl basic tools to recognize and execute these special cases.

```
#]add AnalyticDiffEq ModelingToolkit
using ModelingToolkit, AnalyticDiffEq

# Declare variables and parameters
@parameters t
@derivatives D'~t D2''~t
@variables y(t) x(t)

# Define the Differential Equations: Laguerre, Hermite, and Legendre
eqlag = x*D2(y) + (0 + 1 - x) * D(y) + 6*y ~ 0
eqher = D2(y) - 2*x*D(y) + 5*y ~ 0
eqleg = (1-x^2) * D2(y) - 2*x*D(y) + 4*(4+1)*y ~ 0

# Benchmark against known special case solutions
L6(x) = (1/720)*(x^6 - 36x^5 + 450x^4 - 2400x^3 + 5400x^2 - 4320x + 720)
H5(x) = 120x -160x^3 + 32x^5
P4(x) = (1/8)*(35x^4 - 30x^2 + 3)

@test analyze(eqlag)(.3) ≈ L6(.3)
@test analyze(eqher)(.3) ≈ H5(.3)
@test analyze(eqleg)(.3) ≈ P4(.3)
```

The next steps are to integrate this basic feature into the `solve` interface in DifferentialEquations.jl and
fully leverage the special cases for greater precision, accuracy, and speed.

As this is a proof of concept and a work in progress, `AnalyticDiffEq.jl` can only handle Laguerre, Hermite, and Legendre equations up to order (a hardcoded) order 12, and in the case of Laguerre with `α = 0`.

Extending the use cases to classical orthogonal polynomials will eventually be implemented through [OrthogonalPolynomials.jl](https://github.com/miguelraz/OrthogonalPolynomials.jl) or [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl).
