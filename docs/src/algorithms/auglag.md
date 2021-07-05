# Augmented Lagrangian algorithm

## Description

[Percival.jl](https://github.com/JuliaSmoothOptimizers/Percival.jl) is a pure Julia implementation of the augmented Lagrangian algorithm. Both first and second order versions of the algorithm are available.

## Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using Percival.
```julia
import Percival

alg = AugLag()
options = AugLagOptions()
result = optimize(model, alg, x0, options = options)
```
Percival is an optional dependency of Nonconvex so you need to import it in order to use it.

## Construct an instance

To construct an instance of the Ipopt algorithm, use:
```julia
alg = AugAlg()
```

## Options

The options keyword argument to the `optimize` function shown above must be an instance of the `AugLagOptions` struct when the algorihm is an `AugLag`. To specify options use keyword arguments in the constructor of `AugLagOptions`, e.g:
```julia
options = AugLagOptions(first_order = false, rtol = 1e-4)
```
The most important option is `first_order` which is `true` by default. When `first_order` is `true`, the first order augmented Lagrangian algorithm will be used. And when it is `false`, the second order augmented Lagrangian algorithm will be used. Other arguments include:
- `atol`: absolute tolerance in the subproblem optimizer
- `rtol`: relative tolerance in the subproblem optimizer
- `ctol`: absolute feasibility tolerance
- `max_iter`: maximum number of iterations
- `max_time`: maximum time in seconds
- `max_eval`: maximum number of function evaluations

When using the first order augmented Lagrangian and a block constraint (i.e. a constraint function that returns a vector), the use of reverse-mode AD will only require calling the adjoint operator of the block constraint function in order to compute the gradient of the augmented Lagrangian. This is particularly suitable for constraint functions whose Jacobians are expensive but the adjoint operator is relatively inexpensive.
