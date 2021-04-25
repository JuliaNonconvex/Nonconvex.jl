# Nonconvex

[![Actions Status](https://github.com/mohamed82008/Nonconvex.jl/workflows/CI/badge.svg)](https://github.com/mohamed82008/Nonconvex.jl/actions)
[![codecov](https://codecov.io/gh/mohamed82008/Nonconvex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mohamed82008/Nonconvex.jl)


This package implements and wraps a number of nonconvex constrained optimization algorithms and packages making use of automatic differentiation. The following algorithms are implemented:
- `MMA87`: the original method of moving asymptotes
- `MMA02`: the globally convergent method of moving asymptotes

The following packages are wrapped:
- `IpoptAlg`: a wrapper around Ipopt.jl
- `NLoptAlg`: a wrapper around NLopt.jl
- `AugLag`: a wrapper around Percival.jl which implements the augmented Lagrangian algorithm

The method of moving asymptotes algorithms' were generalized to handle infinite variable bounds. In the augmented Lagrangian algorithm, a block constraint can be handled efficiently by defining a custom adjoint rule for the block constraint using `ChainRulesCore.jl`. This custom adjoint will be picked up by `Nonconvex.jl` when calculating the gradient of the augmented Lagrangian.

# Example

## Load the package

```julia
using Nonconvex
```

## Problem definition

```julia
f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0])
add_ineq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))
```

## Block constraints

```julia
m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0])
block_constr = FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2)
add_ineq_constraint!(m, block_constr)
```

## MMA

```julia
alg = MMA87() # or MMA02()
options = Nonconvex.MMAOptions(
    tol = Nonconvex.Tolerance(kkt = 1e-6, f = 0.0), s_init = 0.1,
)
convcriteria = KKTCriteria()
r = Nonconvex.optimize(
    m, alg, [1.234, 2.345], options = options,
    convcriteria = convcriteria,
)
r.minimum
r.minimizer
```

## NLopt

```julia
alg = NLoptAlg(:LD_MMA)
options = Nonconvex.NLoptOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## Augmented Lagrangian / Percival

```julia
alg = AugLag()
options = Nonconvex.AugLagOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## Mixed equality and inequality constraints

```julia
f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0])
add_eq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))
```

## Ipopt

```julia
alg = IpoptAlg()
options = Nonconvex.IpoptOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## MultistartOptimization

```julia
alg = MsOAlg(TikTak(100))
local_method = NLopt.LN_BOBYQA
options = MsOOptions(use_threads = false)
r = get_starting_point(m, alg, local_method = local_method, options = options)
r.minimum
r.minimizer
```
