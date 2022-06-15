# Optimization result

Each algorithm is free to return a different result type from the `optimize` function. However, all the result types have 2 fields:
- `result.minimum`: stores the minimum objective value reached in the optimization
- `result.minimizer`: stores the optimal decision variables reached during optimization

Some result types store additional information returned by the solver, e.g. the convergence status. Please explore the fields of the `result` output from `optimize` and/or check the documentation of the respective algorithm in the algorithms section of the documentation. If you have further questions, feel free to open issues in the `Nonconvex.jl` repository.
