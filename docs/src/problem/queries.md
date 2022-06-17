# Model queries

There are a number of information you can query about the model after constructing it. These can be useful to check that the model was defined correctly or in the post-processing step after the model has been optimized.

## Number of decision variables

To query the number of decision variables in a model, use:
```julia
NonconvexCore.getnvars(model)
```

## Number of constraints

To query the number of inequality constraints in a model, you can use:
```julia
NonconvexCore.getnineqconstraints(model)
```
A vector-valued constraint will be counted only once.

To query the number of equality constraints, you can use:
```julia
NonconvexCore.getneqconstraints(model)
```

To query the number of semidefinite constraints, you can use:
```julia
NonconvexCore.getnsdconstraints(model)
```

To query the total number of constraints in a model, you can use:
```julia
NonconvexCore.getnconstraints(model)
```
This is the sum of all the previous queries.

## Problem dimension

To get a quick overview of the number of constraints and variables in the model, you can use:
```julia
NonconvexCore.getdim(model)
```
which is  short for:
```julia
(NonconvexCore.getnconstraints(model), NonconvexCore.getnvars(model))
```

## Objective and constraint functions

You can get the objective and constraint functions as regular Julia functions which can be evaluated. To get the objective function, you can use:
```julia
obj = NonconvexCore.getobjective(model)
```

To get a function for all the inequality constraints, you can use:
```julia
ineq = NonconvexCore.getineqconstraints(model)
```

To get a function for all the equality constraints, you can use:
```julia
ineq = NonconvexCore.geteqconstraints(model)
```

To get a function for all the semideifnite constraint functions, you can use:
```julia
ineq = NonconvexCore.getsdconstraints(model)
```

## Initial solution

You can the initial solution using:
```julia
x0 = NonconvexCore.getinit(model)
```

## Variables bounds

You can query the maximum variable bounds for all the variables using:
```julia
NonconvexCore.getmax(model)
```

Similarly, you can query the minimum variable bounds for all the variables using:
```julia
NonconvexCore.getmin(model)
```

## Integrality constraints

To get the vector indiciting which variables are integer, you can use:
```julia
model.integer
```
which will be a `BitVector` (similar to a vector of `Bool`) with `true` corresponding to an integer constraint and `false` corresponding to a continuous variable.
