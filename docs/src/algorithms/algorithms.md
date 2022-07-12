# Algorithms

## Overview of algorithms

A summary of all the algorithms available in `Nonconvex` through different packages is shown in the table below. (scroll right to see more columns)

| Algorithm name | Is meta-algorithm? | Algorithm package | Order | Finite bounds | Infinite bounds | Inequality constraints | Equality constraints | Semidefinite constraints | Integer variables |
| ------- | ----------- | ----- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| Method of moving asymptotes (MMA) | ❌ | `NonconvexMMA.jl` (pure Julia) or `NLopt.jl` | 1 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ |
| Primal dual interior point method | ❌ | `Ipopt.jl` | 1 or 2 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
|  DIviding RECTangles algorithm (DIRECT) | ❌ | `NLopt.jl` | 0 | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ |
| Controlled random search (CRS) | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | 
| Multi-Level Single-Linkage (MLSL) | Limited | `NLopt.jl` | Depends on sub-solver | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| StoGo | ❌ | `NLopt.jl` | 1 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| AGS | ❌ | `NLopt.jl` | 0 | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ |
| Improved Stochastic Ranking Evolution Strategy (ISRES) | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| ESCH | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| COBYLA | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| BOBYQA | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| NEWUOA | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Principal AXIS (PRAXIS) | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Nelder Mead | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Subplex | ❌ | `NLopt.jl` | 0 | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| CCSAQ | ❌ | `NLopt.jl` | 1 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ |
| SLSQP | ❌ | `NLopt.jl` | 1 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| TNewton | ❌ | `NLopt.jl` | 1 | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Shifted limited-memory variable-metric | ❌ | `NLopt.jl` | 1 | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Augmented Lagrangian in `NLopt` | Limited | `NLopt.jl` | Depends on sub-solver | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Augmented Lagrangian in `Percival` | ❌ | `Percival.jl` | 1 or 2 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Multiple trajectory search | ❌ | `NonconvexSearch.jl` | 0 | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Branch and bound for mixed integer nonlinear programming | ❌ | `Juniper.jl` | 1 or 2 | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| Sequential polyhedral outer-approximations for mixed integer nonlinear programming | ❌ | `Pavito.jl` | 1 or 2 | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| Evolutionary centers algorithm (ECA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Differential evolution (DE) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Particle swarm optimization (PSO) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Artificial bee colony (ABC) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Gravitational search algorithm (GSA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Simulated annealing (SA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Whale optimization algorithm (WOA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Machine-coded compact genetic algorithm (MCCGA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Genetic algorithm (GA) | ❌ | `Metaheuristics.jl` | 0 | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Nonlinear optimization with the MADS algorithm (NOMAD) | ❌ | `NOMAD.jl` | 0 | ✅ | ✅ | ✅ | Limited | ❌ | ✅ |
| Topology optimization of binary structures (TOBS) | ❌ | `NonconvexTOBS.jl` | 1 | Binary | ❌ | ✅ | ❌ | ❌ | Binary |
| Hyperband | ✅ | `Hyperopt.jl` | Depends on sub-solver | ✅ | ❌ | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver |
| Random search | ✅ | `Hyperopt.jl` | Depends on sub-solver | ✅ | ❌ | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver |
| Latin hypercube search | ✅ | `Hyperopt.jl` | Depends on sub-solver | ✅ | ❌ | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver |
| Surrogate assisted optimization | ✅ | `NonconvexBayesian.jl` | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver |
| Log barrier method for nonlinear semidefinite constraint handling | ✅ | `NonconvexSemidefinite.jl` | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | Depends on sub-solver | ✅ | Depends on sub-solver |


The following is an explanation of all the columns in the table:
- Algorithm name. This is the name of the algorithm and/or its acronym. Some algorithms have multiple variants implemented in their respective packages. When that's the case, the whole family of algorithms is mentioned only once.
- Is meta-algorithm? Some algorithms are meta-algorithms that call a sub-algorithm to do the optimization after transforming the problem. In this case, a lot of the properties of the meta-algorithm are inherited from the sub-algorithm. So if the sub-algorithm requires gradients or Hessians of functions in the model, the meta-algorithm will also require gradients and Hessians of functions in the model. Fields where the property of the meta-algorithm is inherited from the sub-solver are indicated using the "Depends on sub-solver" entry. Some algorithms in `NLopt` have a "Limited" meta-algorithm status because they can only be used to wrap algorithms from `NLopt`.
- Algorithm package. This is the Julia package that either implements the algorithm or calls it from another programming language. `Nonconvex` wraps all these packages using a consistent API while allowing each algorithm to be customized where possible and have its own set of options.
- Order. This is the order of the algorithm. Zero-order algorithms only require the evaluation of the objective and constraint functions, they don't require any gradients or Hessians of objective and constraint functions. First-order algorithms require both the value and gradients of objective and/or constraint functions. Second-order algorithms require the value, gradients and Hessians of objective and/or constraint functions.
- Finite bounds. This is true if the algorithm supports finite lower and upper bound constraints on the decision variables. One special case is the `TOBS` algorithm which only supports binary decision variables so an entry of "Binary" is used instead of true/false.
- Infinite bounds. This is true if the algorithm supports unbounded decision variables either from below, above or both.
- Inequality constraints. This is true if the algorithm supports nonlinear inequality constraints.
- Equality constraints. This is true if the algorithm supports nonlinear equality constraints. Algorithms that only support linear equality constraints are given an entry of "Limited".
- Semidefinite constraints. This is true if the algorithm supports nonlinear semidefinite constraints.
- Integer variables. This is true if the algorithm supports integer/discrete/binary decision variables, not just continuous. One special case is the `TOBS` algorithm which only supports binary decision variables so an entry of "Binary" is used instead of true/false.

## Wrapper packages

The `JuliaNonconvex` organization hosts a number of packages which wrap other optimization packages in Julia or implement their algorithms. The correct wrapper package is loaded using the `Nonconvex.@load` macro with the algorithm or package name. The following is a summary of all the wrapper packages in the `JuliaNonconvex` organization. To view the documentation of each package, click on the blue docs badge in the last column.

| Package | Description | Tests | Coverage |   Docs   |
| ------- | ----------- | ----- | -------- | -------- |
| [NonconvexMMA.jl](https://github.com/JuliaNonconvex/NonconvexMMA.jl) | Method of moving asymptotes implementation in pure Julia | [![Build Status](https://github.com/JuliaNonconvex/NonconvexMMA.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexMMA.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexMMA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexMMA.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](mma.md) |
| [NonconvexIpopt.jl](https://github.com/JuliaNonconvex/NonconvexIpopt.jl) | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexIpopt.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexIpopt.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexIpopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexIpopt.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](ipopt.md) |
| [NonconvexNLopt.jl](https://github.com/JuliaNonconvex/NonconvexNLopt.jl) | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexNLopt.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexNLopt.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexNLopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexNLopt.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](nlopt.md) |
| [NonconvexPercival.jl](https://github.com/JuliaNonconvex/NonconvexPercival.jl) | [Percival.jl](https://github.com/JuliaSmoothOptimizers/Percival.jl) wrapper (an augmented Lagrangian algorithm implementation) | [![Build Status](https://github.com/JuliaNonconvex/NonconvexPercival.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexPercival.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexPercival.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexPercival.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](auglag.md) |
| [NonconvexJuniper.jl](https://github.com/JuliaNonconvex/NonconvexJuniper.jl) | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexJuniper.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexJuniper.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexJuniper.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexJuniper.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](minlp.md) |
| [NonconvexPavito.jl](https://github.com/JuliaNonconvex/NonconvexPavito.jl) | [Pavito.jl](https://github.com/jump-dev/Pavito.jl) wrapper | [![Build Status](https://github.com/JuliaNonconvex/NonconvexPavito.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexPavito.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexPavito.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexPavito.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](minlp.md) |
| [NonconvexSemidefinite.jl](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl) | Nonlinear semi-definite programming algorithm | [![Build Status](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexSemidefinite.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexSemidefinite.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](sdp.md) |
| [NonconvexMultistart.jl](https://github.com/JuliaNonconvex/NonconvexMultistart.jl) | Multi-start optimization algorithms | [![Build Status](https://github.com/JuliaNonconvex/NonconvexMultistart.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexMultistart.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexMultistart.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexMultistart.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](hyperopt.md) |
| [NonconvexBayesian.jl](https://github.com/JuliaNonconvex/NonconvexBayesian.jl) | Constrained Bayesian optimization implementation | [![Build Status](https://github.com/JuliaNonconvex/NonconvexBayesian.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexBayesian.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexBayesian.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexBayesian.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](surrogate.md) |
| [NonconvexSearch.jl](https://github.com/JuliaNonconvex/NonconvexSearch.jl) | Multi-trajectory and local search methods | [![Build Status](https://github.com/JuliaNonconvex/NonconvexSearch.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexSearch.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexSearch.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexSearch.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](mts.md) |
| [NonconvexTOBS.jl](https://github.com/JuliaNonconvex/NonconvexTOBS.jl) | Binary optimization algorithm called "topology optimization of binary structures" ([TOBS](https://www.sciencedirect.com/science/article/abs/pii/S0168874X17305619?via%3Dihub)) which was originally developed in the context of optimal distribution of material in mechanical components. | [![Build Status](https://github.com/JuliaNonconvex/NonconvexTOBS.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexTOBS.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexTOBS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexTOBS.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](tobs.md) |
| [NonconvexMetaheuristics.jl](https://github.com/JuliaNonconvex/NonconvexMetaheuristics.jl) | Metaheuristic gradient-free optimization algorithms as implemented in [`Metaheuristics.jl`](https://github.com/jmejia8/Metaheuristics.jl). | [![Build Status](https://github.com/JuliaNonconvex/NonconvexMetaheuristics.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexMetaheuristics.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexMetaheuristics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexMetaheuristics.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](metaheuristics.md) |
| [NonconvexNOMAD.jl](https://github.com/JuliaNonconvex/NonconvexNOMAD.jl) | [NOMAD algorithm](https://dl.acm.org/doi/10.1145/1916461.1916468) as wrapped in the [`NOMAD.jl`](https://github.com/bbopt/NOMAD.jl). | [![Build Status](https://github.com/JuliaNonconvex/NonconvexNOMAD.jl/workflows/CI/badge.svg)](https://github.com/JuliaNonconvex/NonconvexNOMAD.jl/actions) | [![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexNOMAD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexNOMAD.jl) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](nomad.md) |

