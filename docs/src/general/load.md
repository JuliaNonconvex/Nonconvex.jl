# Loading algorithm packages with `@load`

`Nonconvex` provides the `@load` macro to install and load the algorithm
packages of the ecosystem on demand. The macro takes no keyword arguments;
it accepts one or more algorithm name symbols:

```julia
using Nonconvex
Nonconvex.@load Ipopt

f(x) = (x[1] - 1.0)^2 + 100.0 * (x[2] - x[1]^2)^2
model = Model(f)
result = optimize(model, IpoptAlg(), [-1.0, 1.0])
```

Several packages can be loaded at once:

```julia
Nonconvex.@load Ipopt NLopt Metaheuristics
```

For each name, the package that implements the corresponding algorithm is
loaded with `using` if it is already installed, and installed with `Pkg.add`
automatically if it is not. The macro must be called with the `Nonconvex.`
qualifier (or imported first with `using Nonconvex: @load`).

## Accepted names

| Names | Package |
| --- | --- |
| `MMA`, `GCMMA`, `MMA87`, `MMA02` | `NonconvexMMA.jl` |
| `Ipopt`, `IpoptAlg` | `NonconvexIpopt.jl` |
| `NLopt`, `NLoptAlg` | `NonconvexNLopt.jl` |
| `Juniper`, `JuniperIpopt`, `JuniperIpoptAlg` | `NonconvexJuniper.jl` |
| `Pavito`, `PavitoIpoptCbc`, `PavitoIpoptCbcAlg` | `NonconvexPavito.jl` |
| `Percival`, `AugLag`, `PercivalAlg` | `NonconvexPercival.jl` |
| `AugLag2` | `NonconvexAugLagLab.jl` |
| `Bayesian`, `BayesOpt`, `BayesOptAlg` | `NonconvexBayesian.jl` |
| `SDP`, `Semidefinite`, `SDPBarrier`, `SDPBarrierAlg` | `NonconvexSemidefinite.jl` |
| `Search`, `MTS`, `LS1`, `RCD`, `MTSAlg`, `LS1Alg`, `RCDAlg` | `NonconvexSearch.jl` |
| `Hyperopt`, `Deflated`, `Multistart`, `HyperoptAlg`, `DeflatedAlg` | `NonconvexMultistart.jl` |
| `TOBS` | `NonconvexTOBS.jl` |
| `Metaheuristics` | `NonconvexMetaheuristics.jl` |
| `NOMAD` | `NonconvexNOMAD.jl` |

Any other name raises an error listing the unsupported algorithm. Each of the
loaded packages exports its algorithm type (for example `IpoptAlg`, `NLoptAlg`,
`MMA87` and `MetaheuristicsECA`); see the [Algorithms](../algorithms/algorithms.md)
section for the documentation of each algorithm package.
