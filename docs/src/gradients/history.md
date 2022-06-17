# Storing history of gradients

Often one may want to store intermediate solutions, function values and gradients for visualisation or post-processing. This is currently not possible with `Nonconvex.jl` as not all solvers support a callback mechanism. To workround this, the `TraceFunction` modifier can be used to store input, output and optionally gradient values
during the optimization:
```julia
F = TraceFunction(f; on_call = false, on_grad = true)
```
`F` can now  be used inplace of `f` in objective and/or constraint functions in a `Nonconvex` model. If the `on_call` keyword argument is set to `true` (default is `true`), the input and output values are stored every time the function `F` is called. If the `on_grad` keyword argument is set to `true` (default is `true`), the input, output and gradient values are stored every time the function `F` is differentiated with either `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote.jl`. The history is stored in `F.trace`. The `TraceFunction` modifier can be compsed with other AD-centric function modifiers in `Nonconvex`, e.g. the `sparsify` or `symbolify` function modifiers.
