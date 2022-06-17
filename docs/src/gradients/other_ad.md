# Using other AD packages

`Nonconvex` uses `Zygote` and `ForwardDiff` by default. There are other AD packages in Julia with different tradeoffs. It is possible to use other AD packages to differentiate specific functions in `Nonconvex` using function modifiers.

[`AbstractDifferentiation.jl`](https://github.com/JuliaDiff/AbstractDifferentiation.jl) is a package that defines a unified API for multiple AD packages. Each AD package has a "backend type" in `AbstractDifferentiation`. You can use any `AbstractDifferentiation`-compatible AD package to differentiate specific functions in `Nonconvex`. The list of `AbstractDifferentiation`-compatible AD packages (other than `Zygote`) are:
- [`FiniteDifferences.jl`](https://github.com/JuliaDiff/FiniteDifferences.jl)
- [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl)
- [`ReverseDiff.jl`](https://github.com/JuliaDiff/ReverseDiff.jl)
- [`Tracker.jl`](https://github.com/FluxML/Tracker.jl)

For more on how to construct a backend struct for each AD package, please refer to the README file of the [`AbstractDifferentiation`](https://github.com/JuliaDiff/AbstractDifferentiation.jl) repository.

In order to use a specific `AbstractDifferentiation`-compatible AD package to differentiate a function `f(x...)` used in a `Nonconvex` objective/constraint, you can use the `abstractdiffy` function modifier from `Nonconvex`:
```julia
F = abstractdiffy(f, backend, x...)
F(x...)
```
where `backend` is an `AbstractDifferentiation` backend struct for the desired AD package, and `x` are all the input arguments to `f`.

The following are common `backend` choices:
- `AbstractDifferentiation.ForwardDiffBackend()` for `ForwardDiff`
- `AbstractDifferentiation.FiniteDifferencesBackend()` for `FiniteDifferences`
- `AbstractDifferentiation.ReverseDiffBackend()` for `ReverseDiff`
- `AbstractDifferentiation.TrackerBackend()` for `Tracker`

Note that in order to define such backend type, one must first load the `AbstractDifferentiation` package as well as the AD package to be used, e.g.:
```
using AbstractDifferentiation, ReverseDiff

backend = AbstractDifferentiation.ReverseDiffBackend()
```

Having defined `F` like this, whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F`, the AD package corresponding to the chosen `backend` will be used instead.

To use [`ForwardDiff`](https://github.com/JuliaDiff/ForwardDiff.jl) as the backend of choice, a shortcut is also available using the `forwarddiffy` function modifier instead of the more general `abstractdiffy`:
```julia
F = forwarddiffy(f, x...)
F(x...)
```
which is short for:
```julia
backend = AbstractDifferentiation.ForwardDiffBackend()
F = abstractdiffy(f, backend, x...)
```
