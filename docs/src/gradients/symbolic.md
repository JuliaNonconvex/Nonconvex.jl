# Symbolic differentiation

## Background

For functions, a tractable symbolic gradient/Jacobian/Hessian may exist. [`Symbolics.jl`](https://github.com/JuliaSymbolics/Symbolics.jl) is a symbolic mathematics package in Julia that can uncover the mathematical expression from Julia functions and then symbolically differentiate the resulting expression. Symbolic simplifications and cancellations can sometimes lead to computational savings compared to algorithmic differentiation. Symbolic differentiation can further exploit the sparsity of the graident, Jacobian and/Hessian if one exists.

In `Nonconvex`, you can enforce the use of `Symbolics` to symbolically differentiate specific functions using the `symbolify` function modifier. In particular, the `Symbolics`-derived gradient/Jacobian/Hessian functions will be used whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate the modified function.

## Symbolifying a function

### First order derivatives

In order to force `Nonconvex` to use `Symbolics` when differentiating a function `f(x...)` once, the `symbolify` function modifier can be used:
```julia
F = symbolify(f, x...; hessian = false, sparse = false, simplify = false)
F(x...)
```
where `x` is some sample input arguments to `f`. `F(x...)` can now be used inplace of `f(x...)` in objectives and/or constraints to be differentiated. Whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F` once, the `Symbolics`-derived gradient/Jacobian will now be used.

The `sparse` keyword argument can be set to `true` (default is `false`) to tell `Symbolics` to return a sparse gradient/Jacobian for the function `F`. The `simplify` keyword argument can be set to `true` (default is `false`) to tell `Symbolics` to simplify the mathematical expressions for the gradient/Jacobian functions.

### Second order derivatives

When `hessian = false` (the default value), only the Jacobian/gradient of `F` will be computed with `Symbolics`. In order to use `Symbolics` to differentiate the function `F` twice, you can set `hessian = true`. Setting `hessian = true` will also work for vector-valued functions `f` or for functions `f` that return multiple, non-vector outputs. The `sparse` and `simplify` keyword arguments work the same way when `hessian` is set to `true`.

## Symbolifying a model

Instead of symbolifying one function at a time, the user can instead symbolify an entire `Nonconvex` model including the objective, all the inequality constraint functions, all the equality constraint functions and all the semidefinite constraint functions.
```julia
sym_model = symbolify(model, hessian = true, simplify = true, sparse = true)
```
where `model` is of type `Model` or `DictModel` and `hessian`, `simplify` and `sparse` have the same intepretation from the function symbolification above. `sym_model` can now be optimized using any of the `Nonconvex` algorithms compatible with the model.

By default, the objective and all the constraint functions will be symbolified. To prevent the symbolification of some component of the model, any of the following keyword arguments can be set to `false` (default is `true`):
- `objective = false`
- `ineq_constraints = false`
- `eq_constraints = false`
- `sd_constraints = false`

Setting the `objective`, `ineq_constraints`, `eq_constraints`, and/or `sd_constraints` keyword arguments to `false` (default is `true`) will prevent the symbolification of the objective, all the inequality constraint functions, all the equality constraint functions, and/or all the semidefinite constraint functions respectively.
