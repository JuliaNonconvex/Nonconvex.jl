# Gradients, Jacobians and Hessians

By default, `Nonconvex` uses:
- The reverse-mode automatic differentiation (AD) package, [`Zygote.jl`](https://github.com/FluxML/Zygote.jl), for computing gradients and Jacobians of functions, and
- The forward-mode AD package, [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl), over  `Zygote.jl` for computing Hessians.

However, one can force `Nonconvex` to use other AD packages or even user-defined gradients and Hessians using special function modifiers. Those special function modifiers customize the behaviour of functions without enforcing the same behaviour on other functions. For instance:
- A specific AD package can be used for one constraint function while the default AD packages are used for other functions in the optimization problem.
- The history of gradients of a specific function can be stored without storing all the gradients of all the functions.
- For functions with a sparse Jacobian or Hessian, the sparsity can be used to speedup the AD using sparse, forward-mode AD for these functions.

In some cases, function modifiers can even be composed on top of each other to create more complex behaviours. 

---
In `Nonconvex`, function modifiers modify the behaviour of a function when differentiated once or twice using either `ForwardDiff` or any [`ChainRules`](https://github.com/JuliaDiff/ChainRules.jl)-compatible AD package, such as `Zygote.jl`. The following features are all implemented in [`NonconvexUtils.jl`](https://github.com/JuliaNonconvex/NonconvexUtils.jl) and re-exported from `Nonconvex`.
---

```@contents
Pages = ["user_defined.md", "other_ad.md", "chainrules_fd.md", "sparse.md", "symbolic.md", "implicit.md", "history.md"]
Depth = 3
```
