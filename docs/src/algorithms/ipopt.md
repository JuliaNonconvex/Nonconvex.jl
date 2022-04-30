# Ipopt

## Description

[Ipopt](https://coin-or.github.io/Ipopt) is a well known interior point optimizer developed and maintained by COIN-OR. The Julia wrapper of Ipopt is [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl). Nonconvex allows the use of Ipopt.jl using the `IpoptAlg` algorithm struct. Ipopt can be used as a second order optimizer using the Hessian of the Lagrangian. Alternatively, an [l-BFGS approximation](https://en.wikipedia.org/wiki/Limited-memory_BFGS) of the Hessian can be used instead turning Ipopt into a first order optimizer tha only requires the gradient of the Lagrangian.

## Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using Ipopt.
```julia
using Nonconvex
Nonconvex.@load Ipopt

alg = IpoptAlg()
options = IpoptOptions()
result = optimize(model, alg, x0, options = options)
```

## Construct an instance

To construct an instance of the Ipopt algorithm, use:
```julia
alg = IpoptAlg()
```

## Options

The options keyword argument to the `optimize` function shown above must be an instance of the `IpoptOptions` struct when the algorihm is an `IpoptAlg`. To specify options use keyword arguments in the constructor of `IpoptOptions`, e.g:
```julia
options = IpoptOptions(first_order = false, tol = 1e-4)
```
There are 2 important and special options:
- `first_order`: `true` by default. When `first_order` is `true`, the first order Ipopt algorithm will be used. And when it is `false`, the second order Ipopt algorithm will be used.
- `linear_constraints`:  `false` by default. When `linear_constraints` is `true`, the Jacobian of the constraints will be computed and sparsified once at the beginning. When it is `false`, dense Jacobians will be computed in every iteration.

All the other options that can be set can be found on the [Ipopt options](https://coin-or.github.io/Ipopt/OPTIONS.html) section of Ipopt's documentation.
