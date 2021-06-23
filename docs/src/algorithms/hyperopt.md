# Multi-start optimization

## Description

[Hyperopt.jl](https://github.com/baggepinnen/Hyperopt.jl) is a Julia library that implements a number of hyperparameter optimization algorithms which can be used to optimize the starting point of the optimization.

## Quick start

Given a model `model` and an initial solution `x0`, the following can be used to optimize the model using Hyperopt.
```julia
import Hyperopt

alg = HyperoptAlg(IpoptAlg())
options = HyperoptOptions(sub_options = IpoptOptions(), sampler = GPSampler())
result = optimize(model, alg, x0, options = options)
```
Hyperopt is an optional dependency of Nonconvex so you need to import it in order to use it. `HyperoptAlg` can wrap any other algorithm in Nonconvex, e.g. `IpoptAlg()`. When the algorithm is a `HyperoptAlg`, the `options` keyword argument must of type `HyperoptOptions`. For more on the options available see below.

## Construct an instance

To construct an instance of the Hyperopt + Ipopt algorithm, use:
```julia
alg = HyperoptAlg(IpoptAlg())
```
`HyperoptAlg` can wrap any other algorithm in Nonconvex, e.g. `NLoptAlg(:LD_MMA)` or `AugLag()`.

## Options

The options keyword argument to the `optimize` function shown above must be an instance of the `HyperoptOptions` struct when the algorihm is a `HyperoptAlg`. To specify options, use keyword arguments in the constructor of `HyperoptOptions`. The `sampler` keyword argument determines the sampling algorithm used to propose new starting points in the multi-start procedure. The `sub_options` keyword argument can be used to pass in the options for the sub-optimizer. There are 2 different ways to pass the sub-options depending on the sampler type.

The `sampler` argument can be of type:
1. `RandomSampler`
2. `LHSampler`
3. `CLHSampler`
4. `GPSampler`
5. `Hyperband`

When optimizing the starting point, the upper and lower bounds on the initial solution must be finite, or finite bounds must be passed in to the `options` constructor. All the options that can be passed to the `HyperoptOptions` constructor are listed below:
```@docs
HyperoptOptions
```

### Sampler choice

#### RandomSampler, LHSampler, CLHSampler and GPSampler

All the sampler constructors are functions defined in Nonconvex wrapping the Hyperopt alternatives to define defaults. For `GPSampler`, `Hyperopt.Min` is always used by default in Nonconvex so you should not pass this argument. All the other arguments that can be passed to the sampler constructor can be found in the [Hyperopt documentation](https://github.com/baggepinnen/Hyperopt.jl#details). Example:
```julia
options = HyperoptOptions(sub_options = IpoptOptions(), sampler = GPSampler())
```

#### Hyperband

The [Hyperband algorithm](https://github.com/baggepinnen/Hyperopt.jl#hyperband) in Hyperopt requires a different way to pass in the sub-options. The Hyperband algorithm tries to optimize the allocation of resources. The `sub_options` argument must be a function with input as the "resources" and output as the sub-solver options. The `Hyperband` constructor accepts 3 arguments:
1. The maximum resources `R`
2. `η` which roughly determines the proportion of trials discarded between each round of successive halving
3. `inner` which specifies an inner sampler of type `RandomSampler`, `LHSampler` or `CLHSampler`.

Example:
```julia
options = HyperoptOptions(
    sub_options = max_iter -> IpoptOptions(max_iter = max_iter), 
    sampler = Hyperband(R=100, η=3, inner=RandomSampler()),
)
```
