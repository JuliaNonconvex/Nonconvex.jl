# `Model` definition

To define an empty model, run:
```julia
model =  Model()
```
To specify an objective function `obj` when creating the model, run:
```julia
model = Model(obj)
```
where `obj` is a function that takes a single vector argument.

## Variable definition

### Add a single variable

To add a new variable to a `Model` with lower and upper bounds `lb` and `ub` respectively, use:
```julia
addvar!(model, lb, ub)
```
The variables added will be stacked on top of each other with a linear integer index. The lower and upper bounds for each variable don't have to be numbers, they can be:
1. Dictionaries
2. Structs
3. `Tuple`s
4. `NamedTuple`s
5. Nested data structures

However, the values input cannot be vectors. A vector input has a different interpretation. See the next section.

The types of the lower and upper bounds of each variable must be the same and this type will be assumed to be the type of the decision variable. Different variables can have different types though. Vectorization and de-vectorization are handled automatically by Nonconvex.

To specify an initial value, use the `init` keyword argument. To add an integer constraint on the variable, use the `integer` keyword argument. For example:
```julia
addvar!(model, 0.0, 10.0, init = 1.0, integer = true)
```
`init` must have the same type as the lower and upper bounds.

### Add multiple variables

To add multiple variables simultaneously, pass in a vector of values for the bounds and optionally for the `init` and `integer` keyword arguments.
```julia
addvar!(model, [0.0, 0.0], [10.0, 10.0], init = [1.0, 1.0], integer = [true, false])
```
The elements of the vector can be:
1. Vectors or arrays in general
2. Dictionaries
3. Structs
4. `Tuple`s
5. `NamedTuple`s
6. Nested data structures

Note that the use of vectors as elements is allowed. Similarly, the types of the lower and upper bounds and the initial values must be the same.

## Objective definition

To specify an objective function after creating the model, use:
```julia
set_objective!(model, obj)
```
where `obj` is a function that takes a single vector argument. The vector input to `obj` will be of the same structure, shape and types as the initial solution, lower bound and upper bound vector.

## Inequality constraint definition

To define an inequality constraint `f(x) <= 0`, where `f` is a Julia function that accepts a single input vector, use:
```julia
add_ineq_constraint!(model, f)
```
The vector input to `f` will be of the same structure, shape and types as the initial solution, lower bound and upper bound vector. The function `f` can return:
1. A number, in which case the constraint will be `f(x) <= 0`
2. A vector or array of numbers, in which case the constraint will be applied element-wise `f(x) .<= 0`.
3. An arbitrary container or data structure, in which case the output will be vectorized first and the constraint will be applied element-wise on the vectorized output.

## Equality constraint definition

To define an inequality constraint `f(x) == 0`, where `f` is a Julia function that accepts a single input vector, use:
```julia
add_eq_constraint!(model, f)
```
The vector input to `f` will be of the same structure, shape and types as the initial solution, lower bound and upper bound vector. The function `f` can return:
1. A number, in which case the constraint will be `f(x) == 0`
2. A vector or array of numbers, in which case the constraint will be applied element-wise `f(x) .== 0`.
3. An arbitrary container or data structure, in which case the output will be vectorized first and the constraint will be applied element-wise on the vectorized output.

## Changing variable bounds

After defining the variables, it is possible to set the minimum and maximum variable bounds to different variables. This is useful for example in iterative procedures where the bounds are updated and the problem is resolved.

To update the entire vector of minimum bounds, you can use:
```julia
setmin!(model, newmin)
```
where `newmin` is a vector of bounds. `-Inf` is allowed in the bounds.

To set a new minimum bound for the `i`th variable only, you can use:
```julia
setmin!(model, i, newmin)
```
instead where `newmin` is a minimum bound of the appropriate type depending on the type of the `i`th variable.

Similarly, to update the entire vector of maximum bounds, you can use:
```julia
setmax!(model, newmax)
```
where `newmax` is a vector of bounds. `Inf` is allowed in the bounds.

To set a new maximum bound for the `i`th variable only, you can use:
```julia
setmax!(model, i, newmax)
```
instead where `newmax` is a maximum bound of the appropriate type depending on the type of the `i`th variable.

## Changing integrality constraints

To constrain a variable to be integer or relax the integrality constraint on the `i`th variable, you can use:
```julia
setinteger!(model, i, integer)
```
where `integer` is `true` to constrain the variable or `false` to relax the constraint.
