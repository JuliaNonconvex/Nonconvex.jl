# Nonconvex

[![Actions Status](https://github.com/JuliaNonconvex/Nonconvex.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/Nonconvex.jl/actions)
[![codecov](https://codecov.io/gh/JuliaNonconvex/Nonconvex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/Nonconvex.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaNonconvex.github.io/Nonconvex.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaNonconvex.github.io/Nonconvex.jl/dev)

`Nonconvex.jl` is an umbrella package over implementations and wrappers of a number of nonconvex constrained optimization algorithms and packages making use of automatic differentiation. Zero, first and second order methods are available. Nonlinear equality and inequality constraints as well as integer constraints are supported. A detailed description of all the algorithms and features available in `Nonconvex` can be found in the [documentation](https://JuliaNonconvex.github.io/Nonconvex.jl/stable).

## The `JuliaNonconvex` organization

The `JuliaNonconvex` organization hosts a number of packages which are available for use in `Nonconvex.jl`. The correct package is loaded using the `Nonconvex.@load` macro with the algorithm or package name. See the [documentation](https://JuliaNonconvex.github.io/Nonconvex.jl/stable) for more details. The following is a summary of all the packages in the `JuliaNonconvex` organization.

| Package | Description | Tests | Coverage |
| ------- | ----------- | ----- | -------- |
| [Nonconvex.jl](https://github.com/mohamed82008/Nonconvex.jl) | Umbrella package for nonconvex optimization | [![Actions Status](https://github.com/JuliaNonconvex/Nonconvex.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/Nonconvex.jl/actions) | [![codecov](https://codecov.io/gh/JuliaNonconvex/Nonconvex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/Nonconvex.jl) |
| [NonconvexCore.jl](https://github.com/JuliaNonconvex/NonconvexCore.jl) | All the interface functions and structs | [![Build Status](https://github.com/JuliaNonconvex/NonconvexCore.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexCore.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexCore.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexCore.jl) |
| [NonconvexMMA.jl](https://github.com/JuliaNonconvex/NonconvexMMA.jl) | Method of moving asymptotes implementation | [![Build Status](https://github.com/JuliaNonconvex/NonconvexMMA.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexMMA.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexMMA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexMMA.jl) |
| [NonconvexIpopt.jl](https://github.com/JuliaNonconvex/NonconvexIpopt.jl) | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexIpopt.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexIpopt.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexIpopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexIpopt.jl) |
| [NonconvexNLopt.jl](https://github.com/JuliaNonconvex/NonconvexNLopt.jl) | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexNLopt.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexNLopt.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexNLopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexNLopt.jl) |
| [NonconvexPercival.jl](https://github.com/JuliaNonconvex/NonconvexPercival.jl) | [Percival.jl](https://github.com/JuliaSmoothOptimizers/Percival.jl) wrapper (an augmented Lagrangian algorithm implementation) | [![Build Status](https://github.com/JuliaNonconvex/NonconvexPercival.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexPercival.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexPercival.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexPercival.jl) |
| [NonconvexJuniper.jl](https://github.com/JuliaNonconvex/NonconvexJuniper.jl) | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexJuniper.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexJuniper.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexJuniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexJuniper.jl) |
| [NonconvexPavito.jl](https://github.com/JuliaNonconvex/NonconvexPavito.jl) | [Pavito.jl](https://github.com/jump-dev/Pavito.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexPavito.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexPavito.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexPavito.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexPavito.jl) |
| [NonconvexSemidefinite.jl](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl) | Nonlinear semi-definite programming algorithm | [![Build Status](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexSemidefinite.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexSemidefinite.jl) |
| [NonconvexMultistart.jl](https://github.com/JuliaNonconvex/NonconvexMultistart.jl) | Multi-start optimization algorithms | [![Build Status](https://github.com/JuliaNonconvex/NonconvexMultistart.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexMultistart.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexMultistart.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexMultistart.jl) |
| [NonconvexBayesian.jl](https://github.com/JuliaNonconvex/NonconvexBayesian.jl) | Constrained Bayesian optimization implementation | [![Build Status](https://github.com/JuliaNonconvex/NonconvexBayesian.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexBayesian.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexBayesian.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexBayesian.jl) |
| [NonconvexSearch.jl](https://github.com/JuliaNonconvex/NonconvexSearch.jl) | Multi-trajectory and local search methods | [![Build Status](https://github.com/JuliaNonconvex/NonconvexSearch.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexSearch.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexSearch.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexSearch.jl) |
| [NonconvexAugLagLab.jl](https://github.com/JuliaNonconvex/NonconvexAugLagLab.jl) | Experimental augmented Lagrangian package | [![Build Status](https://github.com/JuliaNonconvex/NonconvexAugLagLab.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexAugLagLab.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexAugLagLab.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexAugLagLab.jl) |
| [NonconvexUtils.jl](https://github.com/JuliaNonconvex/NonconvexUtils.jl) | Some utility functions for automatic differentiation, history tracing, implicit functions and more. | [![Build Status](https://github.com/JuliaNonconvex/NonconvexUtils.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexUtils.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexUtils.jl) |
| [NonconvexTOBS.jl](https://github.com/JuliaNonconvex/NonconvexTOBS.jl) | Binary optimization algorithm called "topology optimization of binary structures" ([TOBS](https://www.sciencedirect.com/science/article/abs/pii/S0168874X17305619?via%3Dihub)) which was originally developed in the context of optimal distribution of material in mechanical components. | [![Build Status](https://github.com/JuliaNonconvex/NonconvexTOBS.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexTOBS.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexTOBS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexTOBS.jl) |

## Design philosophy

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

## Quick example

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

## How to contribute?

**A beginner?** The easiest way to contribute is to read the documentation, test the package and report issues.

**An impulsive tester?** Improving the test coverage of any package is another great way to contribute to the `JuliaNonconvex` org. Check the coverage report of any of the packages above by clicking the coverage badge. Find the red lines in the report and figure out tests that would cover these lines of code.

**An algorithm head?** There are plenty of optimization algorithms that can be implemented and interfaced in `Nonconvex.jl`. You could be developing the next big nonconvex semidefinite programming algorithm right now! Or the next constraint handling method for evolutionary algorithms!

**A hacker?** Let's figure out how to wrap some optimization package in Julia in the unique, simple and nimble `Nonconvex.jl` style.

**A software designer?** Let's talk about design decisions and how to improve the modularity of the ecosystem.

You can always reach out by opening an issue.
