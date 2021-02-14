"""
    getparent(model::AbstractModel)

Returns the parent model of `model`.
"""
function getparent end

"""
    getobjectiveconstraints(model::AbstractModel)

Returns an instance of `AbstractFunction` that when called with input `x` returns the objective and constraint values of `model` at `x` in a single vector.
"""
function getobjectiveconstraints end

"""
    getobjective(model::AbstractModel)

Returns a function that when called with input `x` returns the objective value of `model` at `x`.
"""
function getobjective end

"""
    getineqconstraint(model::AbstractModel, i::Integer)

Returns a function that when called with input `x` returns the constraint violation value of the `i`^th inequality constraint of `model` at `x`.
"""
function getineqconstraint end

"""
    geteqconstraint(model::AbstractModel, i::Integer)

Returns a function that when called with input `x` returns the constraint violation value of the `i`^th equality constraint of `model` at `x`.
"""
function geteqconstraint end

"""
    getineqconstraints(model::AbstractModel)

Returns a function that when called with input `x` returns the constraint violation values of all the inequality constraints of `model` at `x`.
"""
function getineqconstraints end

"""
    geteqconstraints(model::AbstractModel)

Returns a function that when called with input `x` returns the constraint violation values of all the equality constraints of `model` at `x`.
"""
function geteqconstraints end

"""
    get_objective_multiple(model::AbstractModel)

Returns the factor by which the objective of `model` was scaled. This is 1 by default.
"""
function get_objective_multiple end

"""
    set_objective_multiple!(model::AbstractModel, multiple::Number)

Set the objective's scaling factor of `model` to `multiple`.
"""
function set_objective_multiple! end

"""
    getnconstraints(model::AbstractModel)

Returns the total number of constraints in `model`.
"""
function getnconstraints end

"""
    getnineqconstraints(model::AbstractModel)

Returns the number of inequality constraints in `model`.
"""
function getnineqconstraints end

"""
    getneqconstraints(model::AbstractModel)

Returns the number of equality constraints in `model`.
"""
function getneqconstraints end

"""
    getnvars(model::AbstractModel)

Returns the number of variables in `model`.
"""
function getnvars end

"""
    getmin(model::AbstractModel)

Returns the variables' lower bounds of `model`.
"""
function getmin end

"""
    getmax(model::AbstractModel)

Returns the variables' upper bounds of `model`.
"""
function getmax end

"""
    setmin!(model::AbstractModel, min::Union{Number, AbstractVector})
    setmin!(model::AbstractModel, i::Integer, min::Number)

Sets the variables' lower bounds in `model` to `min`. Or if `i` is given, it sets the lower bound of the `i`^th variable in `model` to `min`.
"""
function setmin! end

"""
    setmax!(model::AbstractModel, max::Union{Number, AbstractVector})
    setmax!(model::AbstractModel, i::Integer, max::Number)

Sets the variables' upper bound in `model` to `max`. Or if `i` is given, it sets the upper bound of the `i`^th variable in `model` to `max`.
"""
function setmax! end

"""
    setbox!(model::AbstractModel, min::Union{Number, AbstractVector}, max::Union{Number, AbstractVector})
    setbox!(model::AbstractModel, i::Integer, min::Number, max::Number)

Sets the variables' lower and upper bounds in `model` to `min` and `max` respectively. Or if `i` is given, it sets the lower and upper bounds of the `i`^th variable in `model` to `min` and `max` respectively.
"""
function setbox! end

"""
    addvar!(model::AbstractModel, lb::Number, ub::Number)

Adds a new variable to `model` with lower and upper bounds `lb` and `ub`.
"""
function addvar! end

"""
    add_ineq_constraint!(model::AbstractModel, f::IneqConstraint)
    add_ineq_constraint!(model::AbstractModel, fs::Vector{<:IneqConstraint})

Adds a new constraint `f` or multiple new constraints `fs` to `model`.
"""
function add_ineq_constraint! end

"""
    init!(model::AbstractModel)

Initializes the model.
"""
function init! end
