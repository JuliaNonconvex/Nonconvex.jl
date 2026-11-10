# The JuliaNonconvex ecosystem

The JuliaNonconvex organization hosts a layered ecosystem of packages built
around a single core:

## [`NonconvexCore.jl`](https://github.com/JuliaNonconvex/NonconvexCore.jl)

The foundation: models, functions and derivatives, the `Workspace`/
`optimize` protocol, solver capability declarations, convergence criteria
and the shared utilities. Every other package in the ecosystem depends on
it and only on it. NonconvexCore 2 is the current major version.

## [`Nonconvex.jl`](https://github.com/JuliaNonconvex/Nonconvex.jl)

The umbrella package (this package). It re-exports NonconvexCore and loads
the solver packages on demand with `Nonconvex.@load` so that a single
`using Nonconvex` gives access to the whole ecosystem.

## Solver packages

Each solver package implements one algorithm family as a
`NonconvexCore.AbstractOptimizer`, declares its exact capabilities, ships
its documentation as a [DocumenterFragments](https://github.com/JuliaDocs/DocumenterFragments.jl)
fragment (integrated into this site) and follows the same quality setup
(TestItems tests, Aqua and JET checks, JuliaFormatter, codecov):

| Package | Algorithm family | Fragment |
|---|---|---|
| `NonconvexIpopt` | interior point (Ipopt) | [Ipopt](ipopt) |
| `NonconvexMMA` | method of moving asymptotes | [MMA](mma) |
| `NonconvexPercival` | augmented Lagrangian (Percival) | [Percival](percival) |
| `NonconvexNLopt` | NLopt algorithm zoo | [NLopt](nlopt) |
| `NonconvexNOMAD` | mesh adaptive direct search | [NOMAD](nomad) |
| `NonconvexMetaheuristics` | evolutionary metaheuristics | [Metaheuristics](metaheuristics) |
| `NonconvexSearch` | multi-trajectory search, RCD | [Search](search) |
| `NonconvexJuniper` | mixed-integer nonlinear programming | [Juniper](juniper) |
| `NonconvexPavito` | mixed-integer nonlinear programming | [Pavito](pavito) |
| `NonconvexTOBS` | topology optimization | [TOBS](tobs) |
| `NonconvexMultistart` | multi-start and hyper-parameter optimization | [Multistart](multistart) |
| `NonconvexSemidefinite` | semidefinite constraints (barrier and re-parameterization) | [Semidefinite](semidefinite) |
| `NonconvexAugLagLab` | augmented Lagrangian laboratory (ALM, PBM, quasi-Newton) | [AugLagLab](auglag) |

## Auxiliary packages

- **`NonconvexBayesian`** — Bayesian optimization (surrogate models and
  acquisition functions); see the [Bayesian](bayesian) fragment.
- **`ReliabilityOptimization`** — reliability analysis and
  reliability-based optimization (FORM/SORM and transformations); see the
  [Reliability](reliability) fragment.

## How the pieces fit together

```text
                  ┌──────────────────────────────┐
                  │        Nonconvex.jl          │   umbrella: @load + reexports
                  └──────────────┬───────────────┘
                                 │ depends on
                  ┌──────────────▼───────────────┐
                  │       NonconvexCore.jl       │   models, derivatives, protocol,
                  │                             │   capabilities, criteria
                  └──────────────┬───────────────┘
          ┌──────────────────────┼────────────────────────┐
          │                      │                        │
   solver packages        auxiliary packages        third-party wrappers
 (Ipopt, MMA, NLopt,   (Bayesian, Reliability,       (TopOpt, ...)
  Percival, NOMAD,       TOBS, ...)
  Metaheuristics,
  Search, Juniper,
  Pavito, Multistart,
  Semidefinite,
  AugLagLab)
```

A model built with NonconvexCore runs unchanged with any solver package:
the solvers implement the same `optimize(model, alg, x0, options...)`
interface, and every algorithm declares its capabilities so that
incompatible models are rejected up front with a descriptive error.

## Documentation

Every package in the ecosystem ships its documentation as a
DocumenterFragments fragment in its `docs/` directory. Each fragment is
built and previewed independently on the package's own pull requests, and
this site assembles all fragments into a single book.
