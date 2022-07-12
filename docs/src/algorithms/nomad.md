# Nonlinear optimization with the MADS (NOMAD) algorithm for continuous and discrete, constrained optimization

## Description

[NOMAD.jl](https://github.com/bbopt/NOMAD.jl) is an optimization package wrapping the C++ implementation of the [NOMAD algorithm](https://dl.acm.org/doi/10.1145/1916461.1916468). `NonconvexNOMAD` allows the use of `NOMAD.jl` using the the `NOMADAlg` struct. `NOMAD.jl` supports continuous and integer decision variables as well as bounds and inequality constraints. Linear equality constraints are also supported when no integer decision variables are in the model.

## Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using `NOMAD`.
```julia
using Nonconvex
Nonconvex.@load NOMAD

alg = NOMADAlg()
options = NOMADOptions()
result = optimize(model, alg, x0, options = options)
```
`NOMAD` is an optional dependency of Nonconvex so you need to load the package to be able to use it.

## Algorithm types

There are 3 different variants of the `NOMADAlg` struct:
- `NOMADAlg(:explicit)`
- `NOMADAlg(:progressive)`
- `NOMADAlg(:custom)`

The explicit algorithm ensures all the constraints are satisfied at all times removing any infeasible point from the population. The progressive algorithm allows infeasible points to be part of the population but enforces feasibility in a progressive manner. The custom variant allows the use of flags on each constraint to declare it as `:explicit` or `:progressive`. For instance, assume `model` is the `Nonconvex` model and `g1` and `g2` are 2 constraint functions.
```julia
add_ineq_constraint!(model, g1, flags = [:explicit])
add_ineq_constraint!(m, g2, flags = [:progressive])
```
The above code declares the first constraint as explicit and the second as progressive. In other words, every point violating the first constraint will be removed from the population but the second constraint will be more progressively enforced.

## Options

The options keyword argument to the `optimize` function shown above must be an instance of the `NOMADOptions` struct when the algorihm is a `NOMADAlg`. To specify options use keyword arguments in the constructor of `NOMADOptions`, e.g:
```julia
options = NOMADOptions()
```
All the options that can be set can be found in the [`NOMAD.jl` documentation](https://bbopt.github.io/NOMAD.jl/stable/nomadProblem/).
