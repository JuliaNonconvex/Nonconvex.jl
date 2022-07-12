# Topology optimization of binary structures (TOBS), a nonlinear binary optimization heuristic 

## Description

The method of topology optimization of binary structures ([TOBS](https://www.sciencedirect.com/science/article/abs/pii/S0168874X17305619?via%3Dihub)) was originally developed in the context of optimal distribution of material in mechanical components. The TOBS algorithm only supports binary decision variables. The TOBS algorithm is a heuristic that relies on the sequential linearization of the objective and constraint functions, progressively enforcing the constraints in the process. The resulting binary linear program can be solved using any mixed integer linear programming (MILP) solver such `Cbc`. This process is repeated iteratively until convergence. This package implements the heuristic for binary nonlinear programming problems.

## Construct an instance

To construct an instance of the `TOBS` algorithm, use:
```julia
alg = TOBSAlg()
```
When optimizing a model using `TOBSAlg`, all the variables are assumed to be binary if their lower and upper bounds are 0 and 1 respectively even if the `isinteger` flag was not used. If there are variables with other bounds' values, the optimization will give an error.

## Example

In this example, the classic topology optimization problem of minimizing the compliance of the structure subject to a volume constraint. Begin by installing and loading the packages required.

```julia
import Nonconvex
Nonconvex.@load TOBS
using Pkg
Pkg.add("TopOpt")
using TopOpt
```

Define the problem and its parameters using [TopOpt.jl](https://github.com/JuliaTopOpt/TopOpt.jl).

```julia
E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
rmin = 6.0 # filter radius
xmin = 0.001 # minimum density
V = 0.5 # maximum volume fraction
p = 3.0 # SIMP penalty

# Define FEA problem
problem_size = (160, 100) # size of rectangular mesh
x0 = fill(1.0, prod(problem_size)) # initial design
problem = PointLoadCantilever(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)
solver = FEASolver(Direct, problem; xmin=xmin)
TopOpt.setpenalty!(solver, p)
cheqfilter = DensityFilter(solver; rmin=rmin) # filter function
comp = TopOpt.Compliance(problem, solver) # compliance function
```

Define the objective and constraint functions.

```julia
obj(x) = comp(cheqfilter(x)) # compliance objective
constr(x) = sum(cheqfilter(x)) / length(x) - V # volume fraction constraint
```

Finally, define the optimization problem using `Nonconvex.jl` and optimize it.

```julia
m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)
options = TOBSOptions()

r = optimize(m, TOBSAlg(), x0; options=options)
r.minimizer
r.minimum
```

The following is a visualization of the optimization history using this example.

![histories](https://user-images.githubusercontent.com/84910559/164938659-797a6a6d-3518-4f7b-a4ff-24b43b822080.png)

![gif](https://user-images.githubusercontent.com/19524993/167059067-f08502a8-c62d-4d62-a2df-e132efc5e25c.gif)

## Options

The following are the options that can be set by passing them to `TOBSOptions`, e.g. `TOBSOptions(movelimit = 0.1)`.
- `movelimit`: the maximum move limit in each iteration as a ratio of the total number of variables. Default value is 0.1, i.e. a maximum of 10% of the variables are allowed to flip value in each iteration.
- `convParam`: the tolerance value. The algrotihm is said to have converged if the moving average of the relative change in the objective value in the last `pastN` iterations is less than `convParam`. Default value is 0.001.
- `pastN`: the number of past iterations used to compute the moving average of the relative change in the objective value. Default value is 20.
- `constrRelax`: the amount of constraint relaxation applied to the linear approximation in each iteration. This is the relative constraint relaxation if the violation is higher than `constrRelax` and the absolute constraint relaxation otherwise. Default value is 0.1.
- `timeLimit`: the time limit (in seconds) of each MILP solve for the linearized sub-problem. Default value is 1.0.
- `optimizer`: the `JuMP` optimizer type used to solve the MILP sub-problem. Default value is `Cbc.Optimizer`.
- `maxiter`: the maximum number of iterations for the algorithm. Default value is 200.
- `timeStable`: a boolean value that when set to `true` switches on the time stability filter of the objective's gradient, discussed in the paper. Default value is `true`.
