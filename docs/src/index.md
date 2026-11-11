# Nonconvex.jl

![JuliaNonconvex logo](https://raw.githubusercontent.com/JuliaNonconvex/Nonconvex.jl/master/docs/src/assets/logo.png)

[Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) is a Julia toolbox for **gradient-based** and **gradient-free** constrained optimization. It assembles the packages of the [JuliaNonconvex](https://github.com/JuliaNonconvex) ecosystem behind one consistent modelling and optimization interface.

## What is in these docs

- **Core** — the modelling and automatic differentiation core provided by [`NonconvexCore.jl`](https://github.com/JuliaNonconvex/NonconvexCore.jl) and re-exported by `Nonconvex`: [`Model`](core/models.md) and `DictModel` for flat and named variables, objectives, constraints, the solver capability system and the [DifferentiationInterface](https://github.com/JuliaDiff/DifferentiationInterface.jl)-based AD machinery. See [Core](core/index.md).
- **Algorithms** — wrappers around well-known optimizers such as [Ipopt](algorithms/ipopt/index.md), [NLopt](algorithms/nlopt/index.md), [MMA](algorithms/mma/index.md), [Percival](algorithms/auglag/index.md), [Juniper](algorithms/juniper/index.md), [Pavito](algorithms/pavito/index.md), [NOMAD](algorithms/nomad/index.md) and [TOBS](algorithms/tobs/index.md), and pure-Julia methods such as [metaheuristics](algorithms/metaheuristics/index.md), [Bayesian optimization](algorithms/bayesian/index.md), [multi-start](algorithms/hyperopt/index.md), [multiple trajectory search](algorithms/mts/index.md), [semidefinite programming](algorithms/sdp/index.md) and [augmented Lagrangian](algorithms/auglaglab/index.md). See the [algorithm overview](algorithms/algorithms.md) for the full capability table.
- **Gradients, Jacobians and Hessians** — reverse- and forward-mode automatic differentiation, user-defined derivatives, [ChainRules](https://github.com/JuliaDiff/ChainRules.jl)-compatible AD backends, sparsity-aware and finite-difference differentiation. See [Gradients, Jacobians and Hessians](gradients/gradients.md).

## Installation

```julia
using Pkg
Pkg.add("Nonconvex")
```

Algorithm packages are installed on demand with the `@load` macro: `Nonconvex.@load Ipopt` installs and loads `NonconvexIpopt.jl`, and likewise for `MMA`, `NLopt`, `Bayesian`, `Metaheuristics`, and the rest of the algorithm packages.

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

- [Core](core/index.md) — models, objectives, constraints and derivatives.
- [Algorithms](algorithms/algorithms.md) — the overview table and the documentation of every algorithm package.
- [Gradients, Jacobians and Hessians](gradients/gradients.md) — how derivatives are computed and customized.
