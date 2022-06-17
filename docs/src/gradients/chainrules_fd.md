# Using ChainRules in ForwardDiff

`ForwardDiff` is a forward-mode AD package that pre-dates `ChainRules`. `ForwardDiff` therefore does not use the `frule`s defined in `ChainRules`. In order to force `ForwardDiff` to use the `frule` defined for a function, one can use the `Nonconvex.NonconvexUtils.@ForwardDiff_frule` macro provided in `Nonconvex`. This is useful in case `ForwardDiff` is used for the entire function but a component of this function has an efficient `frule` defined that you want to take advantage of. To force `ForwardDiff` to use the `frule` defined for a function `f(x::AbstractVector)`, you can use:
```julia
Nonconvex.NonconvexUtils.@ForwardDiff_frule f(x::AbstractVector{<:ForwardDiff.Dual})
```
The signature of the function specifies the method that will be re-directed to use the `frule` from `ChainRules`. Such `frule` therefore needs to be defined for `f` to begin with. `f` with multiple inputs, scalar inputs and other input collection types are also supported.
