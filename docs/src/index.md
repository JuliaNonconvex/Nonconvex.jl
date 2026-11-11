# Nonconvex.jl

![JuliaNonconvex logo](https://raw.githubusercontent.com/JuliaNonconvex/Nonconvex.jl/master/docs/src/assets/logo.png)

[Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) is a Julia toolbox for **gradient-based** and **gradient-free** constrained optimization. It assembles the packages of the [JuliaNonconvex](https://github.com/JuliaNonconvex) ecosystem behind one consistent modelling and optimization interface.

## What is in these docs

- **General** — the modelling and automatic differentiation core provided by [`NonconvexCore.jl`](https://github.com/JuliaNonconvex/NonconvexCore.jl) and re-exported by `Nonconvex`: [`Model`](core/models.md) and `DictModel` for flat and named variables, objectives, constraints, the solver capability system and the [DifferentiationInterface](https://github.com/JuliaDiff/DifferentiationInterface.jl)-based AD machinery. The [Gradients, Jacobians and Hessians](core/gradients/gradients.md) guide covers user-defined derivatives, ChainRules in ForwardDiff and gradient tracing, and [Loading algorithm packages](general/load.md) documents the `@load` macro used to install algorithm packages on demand.
- **Algorithms** — wrappers around well-known optimizers such as [Ipopt](algorithms/ipopt/index.md), [NLopt](algorithms/nlopt/index.md), [MMA](algorithms/mma/index.md), [Percival](algorithms/auglag/index.md), [Juniper](algorithms/juniper/index.md), [Pavito](algorithms/pavito/index.md), [NOMAD](algorithms/nomad/index.md) and [TOBS](algorithms/tobs/index.md), and pure-Julia methods such as [metaheuristics](algorithms/metaheuristics/index.md), [Bayesian optimization](algorithms/bayesian/index.md), [multi-start](algorithms/hyperopt/index.md), [multiple trajectory search](algorithms/mts/index.md), [semidefinite programming](algorithms/sdp/index.md) and [augmented Lagrangian](algorithms/auglaglab/index.md). The full capability table is below.

## Overview of algorithms

A summary of all the algorithms available in `Nonconvex` through different packages is shown in the table below. Scroll right to see more columns and see a description of the columns below the table.

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
## Installation

```julia
using Pkg
Pkg.add("Nonconvex")
```

Algorithm packages are installed on demand with the [`@load`](general/load.md) macro: `Nonconvex.@load Ipopt` installs and loads `NonconvexIpopt.jl`, and likewise for `MMA`, `NLopt`, `Bayesian`, `Metaheuristics`, and the rest of the algorithm packages.

## Quick start

Minimize the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function) using Ipopt:

```julia
using Nonconvex
Nonconvex.@load Ipopt

f(x) = (x[1] - 1.0)^2 + 100.0 * (x[2] - x[1]^2)^2
model = Model(f)
result = optimize(model, IpoptAlg(), [-1.0, 1.0])

result.minimum    # ≈ 0.0
result.minimizer  # ≈ [1.0, 1.0]
```

A `DictModel` works the same way with named variables, and JuMP models can be converted with `DictModel(jump_model)` when JuMP is loaded. Inequality, equality and semidefinite constraints are added with `add_ineq_constraint!`, `add_eq_constraint!` and `add_sd_constraint!` respectively.

## Where to go next

- [General](core/index.md) — models, objectives, constraints and derivatives, plus the [`@load` macro](general/load.md).
- [Algorithms](algorithms/algorithms.md) — the overview table and the documentation of every algorithm package.
- [Gradients, Jacobians and Hessians](core/gradients/gradients.md) — how derivatives are computed and customized.
