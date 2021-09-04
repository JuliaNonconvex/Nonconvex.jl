# Method of moving asymptotes (MMA)

## Description

There are 2 versions of MMA that are available in Nonconvex.jl:
1. The original MMA algorithm from the [1987 paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207).
2. The globally convergent MMA (GCMMA) algorithm from the [2002 paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822).

The MMA algorithm only supports inequality constraints. However, the original algorithm was slightly generalized to handle infinite variable bounds.

## Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using MMA.
```julia
using Nonconvex
Nonconvex.@load MMA

alg = MMA87() # or MMA02()
options = MMAOptions()
result = optimize(model, alg, x0, options = options, convcriteria = KKTCriteria())
```

## Construct an instance

To construct an instance of the original MMA algorithm, use:
```julia
alg = MMA87()
```
or alternatively:
```julia
alg = MMA()
```

To construct an instance of the globally convergent MMA algorithm, use:
```julia
alg = MMA02()
```
or alternatively:
```julia
alg = GCMMA()
```

```@docs
MMA87
MMA02
```

## Options

To specify options for the MMA algorithm, you can construct an instance of `MMAOptions` and use keyword arguments.
```@docs
MMAOptions
```
The `tol` option in MMA can be set to an instance of the `Tolerance` struct:
```@docs
Tolerance
Nonconvex.ConvergenceState
```

## Convergence criteria

There are 4 convergence criteria available for the MMA algorithm:
- `GenericCriteria`
- `KKTCriteria`
- `ScaledKKTCriteria`
- `IpoptCriteria`

```@docs
Nonconvex.ConvergenceCriteria
Nonconvex.GenericCriteria
Nonconvex.KKTCriteria
Nonconvex.ScaledKKTCriteria
Nonconvex.IpoptCriteria
Nonconvex.assess_convergence!
```

To specify the convergence criteria, use:
```julia
converiteria = GenericCriteria()
```
replacing `GenericCriteria()` by `KKTCriteria()`, `ScaledKKTCriteria()` or `IpoptCriteria()`.
