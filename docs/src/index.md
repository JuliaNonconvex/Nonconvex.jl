# Nonconvex.jl Documentation

Nonconvex.jl is a Julia package that implements and wraps a number of constrained nonlinear and mixed integer nonlinear programming solvers. There are 4 unique features of Nonconvex.jl compared to similar packages such as JuMP.jl and NLPModels.jl:

1. Emphasis on a function-based API. Objectives and constraints are normal Julia functions.
2. The use of Zygote.jl for automatic differentiation (AD) of the objective and constraint functions. Specifying analytic gradients is also possible by defining a custom chain rule using ChainRulesCore.jl.
3. The ability to nest algorithms to create more complicated algorithms.
4. The ability to automatically handle structs and different container types in the decision variables by automatically vectorizing and un-vectorizing them in an AD compatible way.

## Installing Nonconvex

To install Nonconvex.jl, open a Julia REPL and type `]` to enter the package mode. Then run:
```julia
add Nonconvex
```

Alternatively, copy and paste the following code to a Julia REPL:
```julia
using Pkg; Pkg.add("Nonconvex")
```

## Loading Nonconvex

To load and start using Nonconvex.jl, run:
```julia
using Nonconvex
```

## Quick start

```julia
using Nonconvex
Nonconvex.@load NLopt

f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

model = Model(f)
addvar!(model, [0.0, 0.0], [10.0, 10.0])
add_ineq_constraint!(model, x -> g(x, 2, 0))
add_ineq_constraint!(model, x -> g(x, -1, 1))

alg = NLoptAlg(:LD_MMA)
options = NLoptOptions()
r = optimize(model, alg, [1.0, 1.0], options = options)
r.minimum # objective value
r.minimzer # decision variables
```

## Table of contents

```@contents
Pages = ["problem/problem.md", "algorithms/algorithms.md", "gradients/gradients.md", "result.md"]
Depth = 3
```
