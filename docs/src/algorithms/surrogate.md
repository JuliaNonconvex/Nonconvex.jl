# Surrogate-assisted Bayesian optimization

## Description

Surrogate-assisted optimization replaces expensive functions in the objecitve and/or constraints by a surrogate. In Nonconvex, a Gaussian process (GP) from [AbstractGPs.jl](https://github.com/JuliaGaussianProcesses/AbstractGPs.jl) is used. A certain amount of "benefit of the doubt" is given to solutions by minimizing:
```julia
μ(x) - η * σ(x)
```
where `μ(x)` and `σ(x)` are the mean and standard deviation of the posterior GP's prediction of the function's value at point `x`.
`η` is a positive number that resembles how much benefit of the doubt we want to give the solution.
A high `η` means more exploration and a low `η` means more exploitation.

Similarly, expensive inequality constraints are replaced by:
```julia
μ(x) - η * σ(x) <= 0
```
giving the solution the benefit of the doubt. And each equality constraint is replaced by 2 inequality constraints as such:
```julia
μ(x) - η * σ(x) <= 0 <= μ(x) + η * σ(x)
```
Once the surrogates are formed, they are solved using a sub-optimizer to get the next query point to update the surrogate model. Prior to the optimization loop, initialization is done using a number of points using a Sobol sequence of points.

## Quick start

```julia
f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

model = Model()
set_objective!(model, f, flags = [:expensive])
addvar!(model, [1e-4, 1e-4], [10.0, 10.0])
add_ineq_constraint!(model, x -> g(x, 2, 0), flags = [:expensive])
add_ineq_constraint!(model, x -> g(x, -1, 1))

alg = BayesOptAlg(IpoptAlg())
options = BayesOptOptions(
    sub_options = IpoptOptions(),
    maxiter = 50, ftol = 1e-4, ctol = 1e-5,
)
r = Nonconvex.optimize(model, alg, [1.234, 2.345], options = options)
```
Note that the `flags` keyword argument was used when defining the objective and constraints and set to `[:expensive]`. This is a hint to Nonconvex to use a surrogate in place of these constraint functions.

## Construct an instance

To construct an instance of the surrogate-assisted optimization algorithm, use:
```julia
alg = BayesOptAlg(subsolver)
```
where `subsolver` is any Nonconvex optimizer to be used to solve the surrogate model.

## Options

The options keyword argument to the `optimize` function shown above must be an instance of the `BayesOptOptions` struct when the algorihm is a `BayesOptAlg`. The following options can be set using keyword arguments when constructing `BayesOptOptions`.
- `sub_options`: options for the sub-optimizer
- `maxiter`: the maximum number of iterations in the Bayesian optimization routine
- `initialize`: `true` by default. If `true`, the GP will be initialized using a Sobol sequence of query points
- `ninit`: number of initialization points
- `ctol`: feasibility tolerance when accepting a solution
- `ftol`: relative tolerance in the function value
- `postoptimize`: `true` by default. If `true`, a local optimization procedure will be used after the Bayesian optimization is completed.
- `kernel`: the GP kernel used. All the kernels from [KernelFunctions.jl](https://github.com/JuliaGaussianProcesses/KernelFunctions.jl) are available.
- `noise`: GP observation noise parameter
- `std_multiple`: `η` in the description of the algorithm above.

## Advanced: manually constructing surrogate functions

Sometimes a function used in the model may need to be replaced by a surrogate but not the entire objective or constraint function. In this case, the surrogate function can be defined explicitly and passed in to the `optimize` function using the keyword argument `surrogates`. A surrogate for the function `f` can be constructed using:
```julia
s1 = Nonconvex.surrogate(f, x0)
```
where `x0` is the initial query point. The output of `s1(x)` will be an interval from [IntervalArithmetic.jl])(https://github.com/JuliaIntervals/IntervalArithmetic.jl) with `lo` and `hi` fields, where `lo = μ(x) - η * σ(x)` and `hi = μ(x) + η σ(x)`. This interval will propagate through the objective function and/or contraint functions outputting an interval or an array of intervals at the end.

To define the objective or constraint functions using the manually contructed surrogates, one needs to return the `lo` field of the output manually at the end of the objective function or inequality constraint function definitions. Equality constraints should also be transformed to a 2-block inequality constraint manually as described above. When manually passing surrogates to the `optimize` function, the `:expensive` flag is redundant and will be ignored.

Example:
```julia
x0 = [1.234, 2.345]
s1 = Nonconvex.surrogate(f, x0)
s2 = Nonconvex.surrogate(x -> [g(x, 2, 0), g(x, -1, 1)], x0)

model = Model()
set_objective!(model, x -> s1(x).lo)
addvar!(model, [1e-4, 1e-4], [10.0, 10.0])
add_ineq_constraint!(model, x -> getproperty.(s2(x), :lo))
alg = BayesOptAlg(IpoptAlg())
options = BayesOptOptions(
    sub_options = IpoptOptions(print_level = 0), maxiter = 50, ctol = 1e-4,
    ninit = 2, initialize = true, postoptimize = false,
)
r = Nonconvex.optimize(model, alg, x0, options = options, surrogates = [s1, s2])
```
