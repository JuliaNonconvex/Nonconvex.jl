# Model definition

There are 2 ways to define a `DictModel`. The direct method is:
```julia
model = DictModel()
```
To pass an objective function while constructing the model, use:
```julia
model = DictModel(obj)
```
where `obj` is a function that takes a single `OrderedDict` argument.

# JuMP model to Nonconvex DictModel

JuMP.jl has an excellent API for defining variables and linear constraints. Using JuMP makes it straightforward to copy a set of linear constraints and variable definitions from a paper. In Nonconvex, you can start with a JuMP model, define variables and constraints using JuMP's API then convert it to a `DictModel`. For example:
```julia
jump_model = JuMP.Model()
@variable jump_model 0 <= x[i=1:3] <= 1
@constraint jump_model sum(x) <= 1
model = DictModel(jump_model)
```
The objective can also be defined either using JuMP or Nonconvex. Once you convert the JuMP model to a Nonconvex model, you can go ahead and define more variables and constraints and/or set the objective in Nonconvex.

# Variable definition

## Add a single variable

Each variable in a `DictModel` has a name which can be a symbol or string. Each named variable can have an arbitrary type, e.g:
1. Vectors or arrays in general
2. Dictionaries
3. Structs
4. `Tuple`s
5. `NamedTuple`s
6. Nested data structures

Vectorization and de-vectorization are handled automatically by Nonconvex.

To add a new named variable to a `DictModel` with a name `:a` and lower and upper bounds `lb` and `ub` respectively, use:
```julia
addvar!(model, :a, lb, ub)
```
Similar to `Model`, optional keyword arguments `init` and `integer` can be set, and the types of the initial value, `lb` and `ub` must be the same.

## Add multiple variables

There is no way to add multiple variables simultaneously to a `DictModel` however a single named variable that's a vector can be added.

# Objective definition

To specify an objective function after creating the model, use:
```julia
set_objective!(model, obj)
```
where `obj` is a function that takes a single `OrderedDict` argument. The `OrderedDict` input to `obj` will be of the same structure, shape and types as the `OrderedDict` initial solution, lower bounds and upper bounds.

# Inequality constraint definition

To define an inequality constraint `f(x) <= 0`, where `f` is a Julia function that accepts a single `OrderedDict` input, use:
```julia
add_ineq_constraint!(model, f)
```
The `OrderedDict` input to `f` will be of the same structure, shape and types as the `OrderedDict` initial solution, lower bounds and upper bounds. The function `f` can return:
1. A number, in which case the constraint will be `f(x) <= 0`
2. A vector or array of numbers, in which case the constraint will be applied element-wise `f(x) .<= 0`.
3. An arbitrary container or data structure, in which case the output will be vectorized first and the constraint will be applied element-wise on the vectorized output.

# Equality constraint definition

To define an inequality constraint `f(x) == 0`, where `f` is a Julia function that accepts a single `OrderedDict` input, use:
```julia
add_eq_constraint!(model, f)
```
The `OrderedDict` input to `f` will be of the same structure, shape and types as the `OrderedDict` initial solution, lower bounds and upper bounds. The function `f` can return:
1. A number, in which case the constraint will be `f(x) == 0`
2. A vector or array of numbers, in which case the constraint will be applied element-wise `f(x) .== 0`.
3. An arbitrary container or data structure, in which case the output will be vectorized first and the constraint will be applied element-wise on the vectorized output.
