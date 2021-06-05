"""
    getparent(model::AbstractModel)

Returns the parent model of `model`.
"""
getparent

"""
    getobjectiveconstraints(model::AbstractModel)

Returns an instance of `AbstractFunction` that when called with input `x` returns the objective and constraint values of `model` at `x` in a single vector.
"""
getobjectiveconstraints

"""
    getobjective(model::AbstractModel)

Returns a that when called with input `x` returns the objective value of `model` at `x`.
"""
getobjective

"""
    getineqconstraint(model::AbstractModel, i::Integer)

Returns a that when called with input `x` returns the constraint violation value of the `i`^th inequality constraint of `model` at `x`.
"""
getineqconstraint

"""
    geteqconstraint(model::AbstractModel, i::Integer)

Returns a that when called with input `x` returns the constraint violation value of the `i`^th equality constraint of `model` at `x`.
"""
geteqconstraint

"""
    getineqconstraints(model::AbstractModel)

Returns a that when called with input `x` returns the constraint violation values of all the inequality constraints of `model` at `x`.
"""
getineqconstraints

"""
    geteqconstraints(model::AbstractModel)

Returns a that when called with input `x` returns the constraint violation values of all the equality constraints of `model` at `x`.
"""
geteqconstraints

"""
    get_objective_multiple(model::AbstractModel)

Returns the factor by which the objective of `model` was scaled. This is 1 by default.
"""
get_objective_multiple

"""
    set_objective_multiple!(model::AbstractModel, multiple::Number)

Set the objective's scaling factor of `model` to `multiple`.
"""
set_objective_multiple!

"""
    getnconstraints(model::AbstractModel)

Returns the total number of constraints in `model`.
"""
getnconstraints

"""
    getnineqconstraints(model::AbstractModel)

Returns the number of inequality constraints in `model`.
"""
getnineqconstraints

"""
    getneqconstraints(model::AbstractModel)

Returns the number of equality constraints in `model`.
"""
getneqconstraints

"""
    getnvars(model::AbstractModel)

Returns the number of variables in `model`.
"""
getnvars

"""
    getmin(model::AbstractModel)

Returns the variables' lower bounds of `model`.
"""
getmin

"""
    getmax(model::AbstractModel)

Returns the variables' upper bounds of `model`.
"""
getmax

"""
    setmin!(model::AbstractModel, min::Union{Number, AbstractVector})
    setmin!(model::AbstractModel, i::Integer, min::Number)

Sets the variables' lower bounds in `model` to `min`. Or if `i` is given, it sets the lower bound of the `i`^th variable in `model` to `min`.
"""
setmin!

"""
    setmax!(model::AbstractModel, max::Union{Number, AbstractVector})
    setmax!(model::AbstractModel, i::Integer, max::Number)

Sets the variables' upper bound in `model` to `max`. Or if `i` is given, it sets the upper bound of the `i`^th variable in `model` to `max`.
"""
setmax!

"""
    setbox!(model::AbstractModel, min::Union{Number, AbstractVector}, max::Union{Number, AbstractVector})
    setbox!(model::AbstractModel, i::Integer, min::Number, max::Number)

Sets the variables' lower and upper bounds in `model` to `min` and `max` respectively. Or if `i` is given, it sets the lower and upper bounds of the `i`^th variable in `model` to `min` and `max` respectively.
"""
setbox!

"""
    addvar!(model::AbstractModel, lb::Number, ub::Number; init = lb, integer = false)

Adds a new variable to `model` with lower and upper bounds `lb` and `ub`, initial value `init` and `integer = true` if the variable is integer.
"""
addvar!

"""
    add_ineq_constraint!(model::AbstractModel, f::IneqConstraint)
    add_ineq_constraint!(model::AbstractModel, fs::Vector{<:IneqConstraint})

Adds a new constraint `f` or multiple new constraints `fs` to `model`.
"""
add_ineq_constraint!

"""
    init!(model::AbstractModel)

Initializes the model.
"""
init!
