# Multi-trajectory search algorithm in pure Julia

## Description

Multiple trajectory search (MTS) is a derivative-free heuristic optimization method presented by [Lin-Yu Tseng and Chun Chen, 2008](https://sci2s.ugr.es/sites/default/files/files/TematicWebSites/EAMHCO/contributionsCEC08/tseng08mts.pdf). 
The `MTS` algorithm is implemented in the `NonconvexSearch.jl` package. This module implements all the optimization methods in the paper.

## Quick start

Using default `MTSOptions()`. `MTS` is used for optimization. 

```julia
using Nonconvex
Nonconvex.@load MTS

alg = MTSAlg()
LS1_options = MTSOptions()
m = Model(f)
lb = [0, 0]
ub = [5, 5]
addvar!(m, lb, ub)
result = optimize(model, alg, x0, options = options)
```

# Random coordinate descent meta-algorithm

## Description

Random coordinate descent (RCD) is a meta-algorithm that optimizes a random
subset of the coordinates in every iteration using a sub-algorithm while
keeping the other coordinates fixed. Breaking the coupling between
coordinates this way can help counter non-convexity. The method is available
in `NonconvexSearch.jl` as `RCDAlg(sub_alg)` where `sub_alg` is any
`AbstractOptimizer` from the Nonconvex ecosystem.

## Quick start

```julia
using Nonconvex
Nonconvex.@load MMA Search

f(x) = sum((x .- 1.0) .^ 2)
m = Model(f)
addvar!(m, fill(-5.0, 4), fill(5.0, 4))

alg = RCDAlg(MMA02())
options = RCDOptions(maxiter = 20, ncoords = 2)
r = optimize(m, alg, fill(0.0, 4), options = options)
```

The options are:

- `maxiter`: the maximum number of outer iterations.
- `fraction`: the fraction of coordinates updated in each iteration. Ignored
  when `ncoords` is set.
- `ncoords`: the number of coordinates updated in each iteration. When
  `nothing`, `fraction` is used.
- `sub_options`: the options passed to the sub-solver. When `nothing`, the
  sub-solver's default options are used.
- `verbose`: whether to print the objective in every iteration.
