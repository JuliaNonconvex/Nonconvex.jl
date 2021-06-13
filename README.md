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
- `JuniperIpoptAlg`: a wrapper around Juniper.jl using Ipopt.jl as a sub-solver
- `PavitoIpoptCbcAlg`: a wrapper around Pavito.jl using Ipopt.jl and Cbc.jl as sub-solvers

The method of moving asymptotes algorithms' were generalized to handle infinite variable bounds. In the augmented Lagrangian algorithm, a block constraint can be handled efficiently by defining a custom adjoint rule for the block constraint using `ChainRulesCore.jl`. This custom adjoint will be picked up by `Nonconvex.jl` when calculating the gradient of the augmented Lagrangian.

# Examples

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
block_constr = x -> [g(x, 2, 0), g(x, -1, 1)]
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
import NLopt

alg = NLoptAlg(:LD_MMA)
options = Nonconvex.NLoptOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## Augmented Lagrangian / Percival

```julia
import Percival

alg = AugLag()
options = Nonconvex.AugLagOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## Equality constraints

You can add equality constraints to the model using:
```julia
add_eq_constraint!(m, f)
```
where `f` is the constraint function to be equal 0.

## Ipopt

```julia
import Ipopt

alg = IpoptAlg()
options = Nonconvex.IpoptOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

## Mixed integer optimization with Juniper and Ipopt

### Juniper

To do mixed integer optimization using Juniper and Ipopt, you can use:
```julia
import Juniper
using Nonconvex

f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0], integer = [false, true])
add_ineq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))

alg = JuniperIpoptAlg()
options = Nonconvex.JuniperIpoptOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer # [0.3327, 1]
```

### Pavito

Or use Pavito with Ipopt and Cbc as sub-solvers using:

```julia
import Pavito
using Nonconvex

f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0], integer = [false, true])
add_ineq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))

alg = PavitoIpoptCbcAlg()
options = Nonconvex.PavitoIpoptCbcOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer # [0.4934, 1.0]
```

## Starting point optimization

### `RandomSampler`, `LHSampler`, `CLHSampler` or `GPSampler`

You can optimize the initial point `x0` using [`Hyperopt.jl`](https://github.com/baggepinnen/Hyperopt.jl):

```julia
import Hyperopt

sampler = RandomSampler()
options = HyperoptOptions(
    sub_options = IpoptOptions(first_order = true),
    sampler = sampler,
)
r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
```

When optimizing the starting point, the upper and lower bounds on the initial solution must be finite. To see all the options that can be set in `HyperoptOptions`, see `?HyperoptOptions`. The sampler can be replaced by `LHSampler()`, `CLHSampler()` or `GPSampler()`. See the documentation of [`Hyperopt.jl`](https://github.com/baggepinnen/Hyperopt.jl) for more details on the options available for each sampler. For `GPSampler`, `Hyperopt.Min` is always used by default in `Nonconvex.jl` so you should not pass this argument.

### `Hyperband`

Alternatively, the hyperband algorithm from `Hyperopt.jl` can be used where the inner sampler can be of type: `RandomSampler`, `LHSampler` or `CLHSampler`:
```julia
sampler = Hyperband(R=100, η=3, inner=RandomSampler())
options = HyperoptOptions(
    sub_options = max_iter -> IpoptOptions(first_order = true, max_iter = max_iter),
    sampler = spl,
)
r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
```
The `sub_options` keyword argument must be a function here that specifies the resources option for the sub-solver used.

## Custom gradient / adjoint

A custom gradient rule for a function should be defined using ChainRulesCore's `rrule`.
For example the following can be used for the function `f` defined above.

```julia
using ChainRulesCore

function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    grad = [0.0, 1 / (2 * sqrt(x[2]))]
    val, Δ -> (NO_FIELDS, Δ * grad)
end
```

You can check it is correct in your tests using [ChainRulesTestUtils.jl](https://github.com/JuliaDiff/ChainRulesTestUtils.jl/).
```julia
using ChainRulesTestUtils
test_rrule(f, [1.2, 3.6])
```

For full details on `rrules` etc see the [ChainRules documentation](https://juliadiff.org/ChainRulesCore.jl/stable/).

## Hack to use other automatic differentiation backends

For specific functions, if you want to use `ForwardDiff` instead of `Zygote`, one way to do this is to define an `rrule` using `ForwardDiff` to compute the gradient or jacobian, e.g:

```julia
using ChainRulesCore, ForwardDiff

function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    grad = ForwardDiff.gradient(f, x)
    val, Δ -> (NO_FIELDS, Δ * grad)
end
```
