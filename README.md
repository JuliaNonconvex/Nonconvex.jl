# Nonconvex

[![Actions Status](https://github.com/mohamed82008/Nonconvex.jl/workflows/CI/badge.svg)](https://github.com/mohamed82008/Nonconvex.jl/actions)
[![codecov](https://codecov.io/gh/mohamed82008/Nonconvex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mohamed82008/Nonconvex.jl)


This package implements and wraps a number of nonconvex constrained optimization algorithms and packages making use of automatic differentiation. The following algorithms are implemented:
- `MMA87`: the original method of moving asymptotes
- `MMA02`: the globally convergent method of moving asymptotes
- `AugLag`: a first-order augmented Lagrangian algorithm

The method of moving asymptotes algorithms' were generalized to handle infinite variable bounds. The augmented Lagrangian implementation may need fine tuning so use with care. In the augmented Lagrangian implementation, a block constraint can be handled efficiently by defining a custom adjoint rule for the block constraint using `ChainRulesCore.jl`. This custom adjoint will be picked up by `Nonconvex.jl` when calculating the gradient of the augmented Lagrangian.

The following packages are wrapped:
- `IpoptAlg`: a wrapper around Ipopt.jl
- `NLoptAlg`: a wrapper around NLopt.jl
- `PercivalAlg`: a wrapper around Percival.jl

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

## Augmented Lagrangian

```julia
alg = AugLag()
options = Nonconvex.AugLagOptions(alg)
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

# Percival

```julia
alg = PercivalAlg()
options = Nonconvex.PercivalOptions()
r = optimize(m, alg, [1.234, 2.345], options = options)
r.minimum
r.minimizer
```

# Hyperopt

You can automatically search a good hyperparameter, using methods in Hyperopt.jl.

Currently support search starting point ```x0```.
```julia
# Automatically search a good starting point (using default option) by add an @search_x0 before optimize call. 
r1 = @search_x0 Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)

# If you prefer customized options
hyperopt_options = X0OptOptions(x0_lb=[0.5, 0.5], x0_rb=[2.8, 2.8],
                                searchspace_size=1000, iters=20, 
                                sampler=Hyperopt.RandomSampler(), 
                                verbose=true,
                                keepall=true)
# Searching hyperparameters using customized options. 
r2 = @hypersearch hyperopt_options, Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
# Equivalent as above. 
r3 = @search_x0 hyperopt_options, Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)

println(r1.minimum)
println(r2.minimum)
println(r3.minimum)
```
