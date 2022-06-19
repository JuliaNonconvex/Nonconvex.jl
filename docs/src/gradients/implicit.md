# Implicit differentiation

## Background

Differentiating implicit functions efficiently using the implicit function theorem has many applications including:
- Nonlinear partial differential equation constrained optimization
- Differentiable optimization layers in deep learning (aka deep declarative networks)
- Differentiable fixed point iteration algorithms for optimal transport (e.g. the Sinkhorn methods)
- Gradient-based bi-level and robust optimization (aka anti-optimization)
- Multi-parameteric programming (aka optimization sensitivity analysis)

For more on implicit differentation, refer to the last part of the [_Understanding automatic differentiation (in Julia)_](https://www.youtube.com/watch?v=UqymrMG-Qi4) video on YouTube and the [_Efficient and modular implicit differentiation_](https://arxiv.org/abs/2105.15183) manuscript for an introduction to the methods implemented here.

## Relationship to [`ImplicitDifferentiation.jl`](https://github.com/gdalle/ImplicitDifferentiation.jl)

[`ImplicitDifferentiation.jl`](https://github.com/gdalle/ImplicitDifferentiation.jl) is an attempt to simplify the implementation in `Nonconvex` making it more lightweight and better documented. For instance, the [documentation of `ImplicitDifferentiation`](https://gdalle.github.io/ImplicitDifferentiation.jl/) presents a number of examples of implicit functions all of which can be defined using `Nonconvex` instead.

## Explicit parameters

There are 4 components to any implicit function:
1. The parameters `p`
2. The variables `x`
3. The residual `f(p, x)` which is used to define `x(p)` as the `x` which satisfies `f(p, x) == 0` for a given value `p`
4. The algorithm used to evaluate `x(p)` satisfying the condition `f(p, x) == 0`

In order to define a differentiable implicit function using `Nonconvex`, you have to specify the "forward" algorithm which finds `x(p)`. For instance, consider the following example:
```julia
using SparseArrays, NLsolve, Zygote, Nonconvex

N = 10
A = spdiagm(0 => fill(10.0, N), 1 => fill(-1.0, N-1), -1 => fill(-1.0, N-1))
p0 = randn(N)

f(p, x) = A * x + 0.1 * x.^2 - p
function forward(p)
  # Solving nonlinear system of equations
  sol = nlsolve(x -> f(p, x), zeros(N), method = :anderson, m = 10)
  # Return the zero found (ignore the second returned value for now)
  return sol.zero, nothing
end
```
`forward` above solves for `x` in the nonlinear system of equations `f(p, x) == 0` given the value of `p`. In this case, the residual function is the same as the function `f(p, x)` used in the forward pass. One can then use the 2 functions `forward` and `f` to define an implicit function using:
```julia
imf = ImplicitFunction(forward, f)
xstar = imf(p0)
```
where `imf(p0)` solves the nonlinear system for `p = p0` and returns the zero `xstar` of the nonlinear system. This function can now be part of any arbitrary Julia function differentiated by Zygote, e.g. it can be part of an objective function in an optimization problem using gradient-based optimization:
```julia
obj(p) = sum(imf(p))
g = Zygote.gradient(obj, p0)[1]
```

In the implicit function's adjoint rule definition, the partial Jacobian `∂f/∂x` is used according to the implicit function theorem. Often this Jacobian or a good approximation of it might be a by-product of the `forward` function. For example when the `forward` function does an optimization using a BFGS-based approximation of the Hessian of the Lagrangian function, the final BFGS approximation can be a good approximation of `∂f/∂x` where the residual `f` is the gradient of the Lagrangian function wrt `x`. In those cases, this Jacobian by-product can be returned as the second argument from `forward` instead of `nothing`.

## Implicit parameters

In some cases, it may be more convenient to avoid having to specify `p` as an explicit argument in `forward` and `f`. The following is also valid to use and will give correct gradients with respect to `p`:
```julia
function obj(p)
  N = length(p)
  f(x) = A * x + 0.1 * x.^2 - p
  function forward()
    # Solving nonlinear system of equations
    sol = nlsolve(f, zeros(N), method = :anderson, m = 10)
    # Return the zero found (ignore the second returned value for now)
    return sol.zero, nothing
  end
  imf = ImplicitFunction(forward, f)
  return sum(imf())
end
g = Zygote.gradient(obj, p0)[1]
```
Notice that `p` was not an explicit argument to `f` or `forward` in the above example and that the implicit function is called using `imf()`. Using some explicit parameters and some implicit parameters is also supported.

## Matrix-free linear solver in the adjoint

In the adjoint definition of implicit functions, a linear system:
```julia
(df/dy) * x = v
```
is solved to find the adjoint vector. To solve the system using a matrix-free iterative solver (GMRES by default) that avoids constructing the Jacobian `df/dy`, you can set the `matrixfree` keyword argument to `true` (default is `false`). When set to `false`, the entrie Jacobian matrix is formed and the linear system is solved using LU factorization.

## Arbitrary data structures

Both `p` and `x` above can be arbitrary data structures, not just arrays of numbers.

## Tolerance

The implicit function theorem assumes that some conditions `f(p, x) == 0` is satisfied. In practice, this will only be approximately satisfied. When this condition is violated, the gradient reported by the implicit function theorem cannot be trusted since its assumption is violated. The maximum tolerance allowed to "accept" the solution `x(p)` and the gradient is given by the keyword argument `tol` (default value is `1e-5`). When the norm of the residual function `f(p, x)` is greater than this tolerance, `NaN`s  are returned for the gradient instead of the value computed via the implicit function theorem. If additionally, the keyword argument `error_on_tol_violation` is set to `true` (default value is `false`), an error is thrown if the norm of the residual exceeds the specified tolerance `tol`.
