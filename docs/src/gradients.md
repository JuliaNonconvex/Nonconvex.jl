# Gradients, Jacobians and Hessians

Nonconvex uses [Zygote.jl](https://github.com/FluxML/Zygote.jl) for automatic differentiation (AD) which in turn uses [ChainRules.jl](https://github.com/JuliaDiff/ChainRules.jl) and [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) for the adjoint rule definitions of different functions.

## Using analytic gradients

To use an analytically derived gradient function `analytic_gradient` for the function `f`, you need to define an adjoint rule for the function `f` using `ChainRulesCore.jl` as such:
```julia
function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    grad = analytic_gradient(f, x)
    val, Δ -> (NoTangent(), Δ * grad)
end
```
It's always good practice to check that the rule deifned is correct using [ChainRulesTestUtils.jl](https://github.com/JuliaDiff/ChainRulesTestUtils.jl/).
```julia
using ChainRulesTestUtils
test_rrule(f, [1.2, 3.6])
```
For full details on rrules etc see the [ChainRules documentation](https://juliadiff.org/ChainRulesCore.jl/stable/).

## Using analytic Jacobians

To use an analytically derived jacobian function `analytic_jacobian` for the function `f`, you need to define an adjoint rule for the function `f` using `ChainRulesCore.jl` as such:
```julia
function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    jac = analytic_jacobian(f, x)
    val, Δ -> (NoTangent(), jac' * Δ)
end
```
It's always good practice to check that the rule deifned is correct using [ChainRulesTestUtils.jl](https://github.com/JuliaDiff/ChainRulesTestUtils.jl/).
```julia
using ChainRulesTestUtils
test_rrule(f, [1.2, 3.6])
```
For full details on rrules etc see the [ChainRules documentation](https://juliadiff.org/ChainRulesCore.jl/stable/).

## Using analytic Hessians

There is currently no way to use analytic Hessianas in Nonconvex.jl. If this is a feature you need, please an open an issue and it may be added.

## Using other automatic differentiation backends

For specific functions if you want to use [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) instead of Zygote to compute the gradient or Jacobian, you can define an `rrule` that uses ForwardDiff to compute the gradient or jacobian, e.g:
```julia
using ChainRulesCore, ForwardDiff

function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    grad = ForwardDiff.gradient(f, x)
    val, Δ -> (NoTangent(), Δ * grad)
end
```
or
```julia
using ChainRulesCore, ForwardDiff

function ChainRulesCore.rrule(::typeof(f), x::AbstractVector)
    val = f(x)
    jac = ForwardDiff.jacobian(f, x)
    val, Δ -> (NoTangent(), jac' * Δ)
end
```

## AD limitations in Nonconvex

- Sparse Jacobians are not yet supported
- Analytic Hessian functions are not possible to use
