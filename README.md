# Nonconvex

[![Actions Status](https://github.com/mohamed82008/Nonconvex.jl/workflows/CI/badge.svg)](https://github.com/mohamed82008/Nonconvex.jl/actions)
[![codecov](https://codecov.io/gh/mohamed82008/Nonconvex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mohamed82008/Nonconvex.jl)
[![Documentation](https://img.shields.io/badge/doc-latest-blue.svg)](https://mohamed82008.github.io/Nonconvex.jl/dev)


This package implements and wraps a number of nonconvex constrained optimization algorithms and packages making use of automatic differentiation. The following algorithms are implemented:
- `MMA87`: the original method of moving asymptotes
- `MMA02`: the globally convergent method of moving asymptotes
- Surrogate-assisted optimization using Gaussian processes

The following packages are wrapped:
- `IpoptAlg`: a wrapper around [`Ipopt.jl`](https://github.com/jump-dev/Ipopt.jl)
- `NLoptAlg`: a wrapper around [`NLopt.jl`](https://github.com/JuliaOpt/NLopt.jl)
- `AugLag`: a wrapper around [`Percival.jl`](https://github.com/JuliaSmoothOptimizers/Percival.jl) which implements the augmented Lagrangian algorithm
- `JuniperIpoptAlg`: a wrapper around [`Juniper.jl`](https://github.com/lanl-ansi/Juniper.jl) using [`Ipopt.jl`](https://github.com/jump-dev/Ipopt.jl) as a sub-solver
- `PavitoIpoptCbcAlg`: a wrapper around [`Pavito.jl`](https://github.com/jump-dev/Pavito.jl) using [`Ipopt.jl`](https://github.com/jump-dev/Ipopt.jl) and [`Cbc.jl`](https://github.com/jump-dev/Cbc.jl) as sub-solvers
- `HyperoptAlg`: a wrapper around [`Hyperopt.jl`](https://github.com/baggepinnen/Hyperopt.jl) with your choice of sub-solver from any of the other solvers in this package.

The method of moving asymptotes algorithms' were generalized to handle infinite variable bounds. In the augmented Lagrangian algorithm, a block constraint can be handled efficiently by defining a custom adjoint rule for the block constraint using [`ChainRulesCore.jl`](https://github.com/JuliaDiff/ChainRulesCore.jl). This custom adjoint will be picked up by `Nonconvex.jl` when calculating the gradient of the augmented Lagrangian.

The documentation of all the algorithms and features available in Nonconvex can be found [here](docs/src/index.md).

```@contents
Pages = ["docs/src/problem.md", "docs/src/algorithms/algorithms.md", "docs/src/gradients.md"]
Depth = 1
```
