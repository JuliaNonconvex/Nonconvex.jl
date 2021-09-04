# Mixed integer nonlinear programming (MINLP)

## Description

There are 2 MINLP solvers available in Nonconvex:
1. [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl) with [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) as a sub-solver.
2. [Pavito.jl](https://github.com/jump-dev/Pavito.jl) with [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) and [Cbc.jl](https://github.com/jump-dev/Cbc.jl) as sub-solvers.

These rely on local nonlinear programming solvers and a branch and bound procedure to find a locally optimal solution that satisfies the integerality constraints.

## Juniper + Ipopt

### Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using Juniper and Ipopt.
```julia
using Nonconvex
Nonconvex.@load Juniper

alg = JuniperIpoptAlg()
options = JuniperIpoptOptions()
result = optimize(model, alg, x0, options = options)
```
Juniper is an optional dependency of Nonconvex, so you need to load it in order to use it. Note that the integer constraints must be specified when defining variables. See the [problem definition](../problem.md) documentation for more details.

### Construct an instance

To construct an instance of the Juniper + Ipopt algorithm, use:
```julia
alg = JuniperIpoptAlg()
```

### Options

The options keyword argument to the `optimize` function shown above must be an instance of the `JuniperIpoptOptions` struct when the algorihm is a `JuniperIpoptAlg`. To specify options use, keyword arguments in the constructor of `JuniperIpoptOptions`, e.g:
```julia
options = JuniperIpoptOptions(first_order = false, linear_constraints = true, subsolver_options = IpoptOptions(), atol = 1e-4)
```
There are 3 important and special options you can pass to the optimizer:
- `first_order`: `true` by default. When `first_order` is `true`, the first order Ipopt algorithm will be used. And when it is `false`, the second order Ipopt algorithm will be used.
- `linear_constraints`: `false` by default. When `linear_constraints` is `true`, the Jacobian of the constraints will be computed and sparsified once at the beginning. When it is `false`, dense Jacobians will be computed in every iteration.
- `subsolver_options`: an instance of `IpoptOptions` to be used in the Ipopt sub-solver.

All the other options to Juniper can be found in the [Juniper documentation](https://lanl-ansi.github.io/Juniper.jl/stable/options/).

## Pavito + Ipopt + Cbc

### Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using Juniper and Ipopt.
```julia
import Pavito

alg = PavitoIpoptCbcAlg()
options = PavitoIpoptCbcOptions()
result = optimize(model, alg, x0, options = options)
```
Pavito is an optional dependency of Nonconvex, so you need to load it in order to use it. Note that the integer constraints must be specified when defining variables. See the [problem definition](../problem.md) documentation for more details.

### Construct an instance

To construct an instance of the Pavito + Ipopt + Cbc algorithm, use:
```julia
alg = PavitoIpoptCbcAlg()
```

### Options

The options keyword argument to the `optimize` function shown above must be an instance of `PavitoIpoptCbcOptions` struct when the algorithm is a `PavitoIpoptCbcAlg`. To specify options, use keyword arguments in the constructor of `JuniperIpoptOptions` or `PavitoIpoptCbcOptions`, e.g:
```julia
options = PavitoIpoptCbcOptions(first_order = false, subsolver_options = IpoptOptions(), timeout = 120.0)
```
There are 2 important and special options you can pass to the optimizer:
- `first_order`: `true` by default. When `first_order` is `true`, the first order Ipopt algorithm will be used. And when it is `false`, the second order Ipopt algorithm will be used.
- `subsolver_options`: an instance of `IpoptOptions` to be used in the Ipopt sub-solver.

All the other options to Pavito can be found in the [Pavito documentation](https://github.com/jump-dev/Pavito.jl#pavito-solver-options).
