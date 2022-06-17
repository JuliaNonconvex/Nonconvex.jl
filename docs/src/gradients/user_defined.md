# User-defined gradient, Jacobian or Hessian

## Gradients and Jacobians

To use a user-defined gradient/Jacobian function `g(x)` for a function `f(x)`, you can use the `CustomGradFunction` modifier:
```julia
F = CustomGradFunction(f, g)
F(x)
```
`F` can be then used in place of `f` as an objective function, as a constraint function or as part of any such function. When `f` is scalar-valued, `g` is expected to return a gradient vector. When `f` is vector-valued, `g` is expected to return a Jacobian matrix. Whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F`, the custom gradient/Jacobian function `g` will be used.

## Hessian or Hessian vector product

For second-order optimization algorithms, a user-defined Hessian function `h(x)` can be used for any scalar-valued function `f(x)` with gradient `g(x)`. To use a user-defined Hessian function `h(x)`, you can use the `CustomHessianFunction` modifier:
```julia
F = CustomHessianFunction(f, g, h)
F(x)
```
`F` can be then used in place of `f` as an objective function, as a constraint function or as part of any such function. `f` is expected to return a scalar, `g` is expected to return a gradient vector and `h` is expected to return a symmetric Hessian matrix. Whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F`, the custom gradient and Hessian functions will be used.

Instead of a Hessian function `h`, alternatively a Hessian-vector product operator `h(x, v)` can be used, which multiplies the Hessian of `f` at `x` by the vector `v`. To use a Hessian-vector product operator `hvp(x, v)` instead of computing the full Hessian, you can pass the `hvp` function as the third argument to `CustomHessianFunction` and set the `hvp` keyword argument to `true`:
```julia
F = CustomHessianFunction(f, g, hvp; hvp = true)
F(x)
```
