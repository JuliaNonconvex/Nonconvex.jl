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
