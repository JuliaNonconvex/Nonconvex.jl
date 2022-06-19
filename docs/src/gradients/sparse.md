# Sparse Jacobian or Hessian

## Background

For functions with a sparse Jacobian or Hessian, it can sometimes be useful to exploit such sparsity to speedup the computation of the Jacobian. This can be done using the [`SparseDiffTools.jl`](https://github.com/JuliaDiff/SparseDiffTools.jl) package.

`SparseDiffTools` can compute multiple columns of the Jacobian matrix of a vector-valued function `y = f(x)` simulatenously using a single Jacobian-vector product operation. Such columns corresponding to a subset of the input variables, e.g. `(x[1], x[3])`, however need not overlap in the output variables they influence. For instance, assume
- `y[1]` and `y[2]` are a function of `x[1]` and `x[2]` only, and
- `y[3]` and `y[4]` are a function of `x[3]` and `x[4]` only.

The Jacobian `dy/dx` will therefore have a block diagonal structure. Additionally, since `x[1]` and `x[3]` do not affect the same output variables, their corresponding columns in the block-diagonal Jacobian can be computed simulatenously using a single Jacobian-vector block. The same thing for `x[2]` and `x[4]`. Finding such subsets of input variables such that no 2 input variables in a subset affect the same output is done using `SparseDiffTools`. In the diagonal Jacobian case, all the input variables do not overlap so all the columns of the Jacobian can be obtained using a single Jacobian-vector product.

The problem of finding the optimal splitting of input variables to require the least number of Jacobian-vector products when computing the full Jacobian can be formulated as a [graph coloring](https://en.wikipedia.org/wiki/Graph_coloring) problem in computer science, which is an [NP-hard](https://en.wikipedia.org/wiki/NP-hardness) problem. `SparseDiffTools` uses a tractable heuristic to find reasonable splittings for different Jacobian or Hessian sparsity patterns. The sparsity pattern of the Jacobian or Hessian can either be user-provided or it will be automatically uncovered using [`Symbolics.jl`](https://github.com/JuliaSymbolics/Symbolics.jl).

In `Nonconvex`, you can enforce the use of `SparseDiffTools` for specific functions using the `sparsify` function modifier. In particular, the `rrule` and `frule` of the modified function will be using `SparseDiffTools` to find the full Jacobian first and then doing either a Jacobian-vector product in the `frule` or a vector-Jacobian product in the `rrule`. Such `frule` will also be used by `ForwardDiff` if used to differentiate the modified function. For more on `frule`s, `rrule`s, Jacobian-vector products and vector-Jacobian products, refer to the following video on [Understanding autmoatic differentiation (in Julia)](https://www.youtube.com/watch?v=UqymrMG-Qi4).

## Sparsifying a function

### First order derivatives

In order to force `Nonconvex` to use `SparseDiffTools` when differentiating a function `f(x...)` once, the `sparsify` function modifier can be used:
```julia
F = sparsify(f, x...; hessian = false)
F(x...)
```
where `x` is some sample input arguments to `f`. `F(x...)` can now be used inplace of `f(x...)` in objectives and/or constraints to be differentiated. Whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F` once, `SparseDiffTools` will now be used.

### Second order derivatives

When `hessian = false` (the default value), only the Jacobian/gradient of `F` will be treated as sparse. In order to use `SparseDiffTools` to compute sparse second order derivatives as well, you can set `hessian = true`. This is recommended for scalar-valued functions with sparse Hessian matrices. Setting `hessian = true` will also work for vector-valued functions `f` or for functions `f` that return multiple, non-vector outputs. The sparsity of the third order derivative tensor will be used to compute the third order tensor efficiently.

### User-defined sparsity patterns

Using `sparsify` as shown above will make use of `Symbolics` to uncover the sparsity of the Jacobian and Hessian matrices of `f`. In some cases, the function `f` may not be `Symbolics`-compatible or it may have a known sparsity pattern. The user can therefore use the `jac_pattern` or `hess_pattern` keyword arguments to set the pattern manually.

The `jac_pattern` is expected to be a `SparseMatrixCSC` of element type `Bool` with `true` where there is a structural non-zero, and `false` where there is a structural zero in the Jacobian matrix. The size of `jac_pattern` should be `noutputs x ninputs` where `noutputs` is the number of outputs of `f` and `ninputs` is the number of inputs to `f`. When the inputs and/or outputs are multiple and/or non-vector, they are assumed to be flattened to vectors and `noutputs`/`ninputs` is the length of the flat vector.

Passing the Hessian sparsity pattern is also possible using the `hess_pattern` keyword argument. For scalar-valued functions, `hess_pattern` should have size `ninputs x ninputs` where `ninputs` is the number of input variables in the flattened input arguments. For vector-valued functions, the sparsity pattern will be the sparsity pattern of the Jacobian of the linearized Jacobian of `f`. Assume `f` takes a single vector input `x` and returns a single output vector `y`. The sparsity pattern will be that of `d(vec(dy/dx))/dx`. `hess_pattern` should therefore have size `(noutputs * ninputs) x ninputs`. For example, assume `y` is a vector of length 2 and `x` is a vector of length 3. The Jacobian `dy/dx` will be a matrix of size `2 x 3`. `vec(dy/dx)` will be a vector of length 6. `d(vec(dy/dx))/dx` will be a matrix of size `6 x 3`. `hess_pattern` should therefore be a `SparseMatrixCSC` with element type `Bool` and size `6 x 3`.

For general functions `f` with multiple or non-vector inputs or outputs, `noutputs` and `ninputs` are the lengths of the flattened outputs and inputs respectively.

## Sparsifying a model

Instead of sparsifying one function at a time, the user can instead sparsify an entire `Nonconvex` model including the objective, all the inequality constraint functions, all the equality constraint functions and all the semidefinite constraint functions.
```julia
sp_model = sparsify(model, hessian = true)
```
where `model` is of type `Model` or `DictModel` and `hessian` has the same intepretation from the function sparisfication above. `sp_model` can now be optimized using any of the `Nonconvex` algorithms compatible with the model.

By default, the objective and all the constraint functions will be sparsified. To prevent the sparisfication of some component of the model, any of the following keyword arguments can be set to `false` (default is `true`):
- `objective = false`
- `ineq_constraints = false`
- `eq_constraints = false`
- `sd_constraints = false`

Setting the `objective`, `ineq_constraints`, `eq_constraints`, and/or `sd_constraints` keyword arguments to `false` (default is `true`) will prevent the sparisification of the objective, all the inequality constraint functions, all the equality constraint functions, and/or all the semidefinite constraint functions respectively.
