# Interior point method using `Ipopt.jl`

## Description

[Ipopt](https://coin-or.github.io/Ipopt) is a well known interior point optimizer developed and maintained by COIN-OR. The Julia wrapper of Ipopt is [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl). `Ipopt.jl` is wrapped in `NonconvexIpopt.jl`. `NonconvexIpopt` allows the use of `Ipopt.jl` using the `IpoptAlg` algorithm struct. `IpoptAlg` can be used as a second order optimizer computing the Hessian of the Lagrangian in every iteration. Alternatively, an [l-BFGS approximation](https://en.wikipedia.org/wiki/Limited-memory_BFGS) of the Hessian can be used instead turning `IpoptAlg` into a first order optimizer tha only requires the gradient of the Lagrangian.

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
options = IpoptOptions(first_order = false, tol = 1e-4, sparse = false, symbolic = false)
```
There are 4 important and special options:
- `first_order`: `true` by default. When `first_order` is `true`, the first order Ipopt algorithm will be used. And when it is `false`, the second order Ipopt algorithm will be used.
- `symbolic`: `false` by default. When `symbolic` is set to `true`, the gradients, Jacobians and Hessians of the objective, constraint and Lagrangian functions will be calculated using symbolic differentiation from [`Symbolics.jl`](https://github.com/JuliaSymbolics/Symbolics.jl). This is the same approach used by `symbolify` which is described in the [symbolic differentiation section](../gradients/symbolic.md) in the documentation.
- `sparse`: `false` by default. When `sparse` is set to `true`, the gradients, Jacobians and Hessians of the objective, constraint and Lagrangian functions will be treated as sparse vectors/matrices. When combined with `symbolic = true`, the output of symbolic differentiation will be a sparse vector/matrix, akin to setting `sparse = true` in the `symbolify` function discussed in [symbolic differentiation section](../gradients/symbolic.md) in the documentation. When used alone with `symbolic = false`, [`SparseDiffTools.jl`](https://github.com/JuliaDiff/SparseDiffTools.jl) is used instead for the differentiation and `Symbolics` is only used to get the sparsity pattern, much like how `sparsify` works. For more details on `sparsify` and the way `SparseDiffTools` works, see the [sparsity section](../gradients/sparse.md) in the documentation is used instead.
- `linear_constraints`:  `false` by default. When `linear_constraints` is `true`, the Jacobian of the constraints will be computed and sparsified once at the beginning. When it is `false`, dense Jacobians will be computed in every iteration.

> Note that there is no need to use `sparsify` or `symbolify` on the model or functions before optimizing it with an `IpoptAlg`. Setting the `sparse` and `symbolic` options above are enough to trigger the symbolic differentiation and/or sparsity exploitation.

All the other options that can be set can be found on the [Ipopt options](https://coin-or.github.io/Ipopt/OPTIONS.html) section of Ipopt's documentation.
