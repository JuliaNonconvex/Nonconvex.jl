# Nonconvex

This package implements and wraps a number of nonconvex constrained optimization algorithms and packages making use of automatic differentiation. The following algorithms are implemented:
- `MMA87`: the original method of moving asymptotes
- `MMA02`: the globally convergent method of moving asymptotes
- `AugLag`: a first-order augmented Lagrangian algorithm

The method of moving asymptotes algorithms' were generalized to handle infinite variable bounds. The augmented Lagrangian implementation may need fine tuning so use with care. In the augmented Lagrangian implementation, a block constraint can be handled efficiently by defining a custom adjoint rule for the block constraint using `ChainRulesCore.jl`. This custom adjoint will be picked up by `Nonconvex.jl` when calculating the gradient of the augmented Lagrangian.

The following packages are wrapped:
- `IpoptAlg`: a wrapping around Ipopt.jl

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

## Augmented Lagrangian

```julia
alg = AugLag()
options = Nonconvex.AugLagOptions(alg)
r = optimize(m, alg, x0, options = options)
```

## Ipopt

```julia
alg = IpoptAlg()
options = Nonconvex.IpoptOptions()
r = optimize(m, alg, x0, options = options)
```
