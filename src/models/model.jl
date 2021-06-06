abstract type AbstractModel end

"""
```
struct Model <: AbstractModel
    objective::Objective
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::AbstractVector
    box_max::AbstractVector
    init::AbstractVector
    integer::AbstractVector
end
```

The `Model` structs stores information about the nonlinear optimization problem.
- `objective`: the objective function of the problem of type [`Objective`](@ref).
- `eq_constraints`: the equality constraints of the problem of type [`VectorOfFunctions`](@ref). Each function in `ineq_constraints` can be an instance of [`IneqConstraint`](@ref) or [`AbstractFunction`](@ref). If the function is not an `IneqConstraint`, its right-hand-side bound is assumed to be 0.
- `ineq_constraints`: the inequality constraints of the problem of type [`VectorOfFunctions`](@ref). Each function in `ineq_constraints` can be an instance of [`IneqConstraint`](@ref) or [`AbstractFunction`](@ref). If the function is not an `IneqConstraint`, its right-hand-side bound is assumed to be 0.
"""
mutable struct Model{Tv <: AbstractVector} <: AbstractModel
    objective::Union{Nothing, Objective}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::Tv
    box_max::Tv
    init::Tv
    integer::BitVector
end

"""
    Model()
    Model(f)

Constructs an empty model or a model with objective function `f`. The decision variables are assumed to be of type `Vector{Any}`.
"""
Model(f::Union{Nothing, Function} = nothing) = Model(Any, f)

"""
    Model(::Type{T}, f::Union{Nothing, Function}) where {T}

Constructs an empty model with objective function `f` and decision variable value type `Vector{T}`. `f` can be an instance of `Base.Function` but must return a number,  or it can be an intance of [`Objective`](@ref). `f` can also be `nothing` in which case, the objective function can be defined later using [`set_objective!`](@ref).
"""
Model(::Type{T}, f::Function) where {T} = Model(T, Objective(f))
function Model(::Type{T}, obj::Union{Nothing, Objective}) where {T}
    return Model(obj, VectorOfFunctions(EqConstraint[]), VectorOfFunctions(IneqConstraint[]), T[], T[], T[], falses(0))
end

getobjective(m::AbstractModel) = m.objective

getineqconstraints(m::AbstractModel) = m.ineq_constraints

geteqconstraints(m::AbstractModel) = m.eq_constraints

function getobjectiveconstraints(m::AbstractModel)
    return VectorOfFunctions((getobjective(m), getineqconstraints(m), geteqconstraints(m)))
end

getineqconstraint(m::AbstractModel, i::Integer) = getineqconstraints(m).fs[i]

geteqconstraint(m::AbstractModel, i::Integer) = geteqconstraints(m).fs[i]

getnineqconstraints(m::AbstractModel) = length(getineqconstraints(m))

getneqconstraints(m::AbstractModel) = length(geteqconstraints(m))

getnconstraints(m::AbstractModel) = getnineqconstraints(m) + getneqconstraints(m)

getnvars(m::AbstractModel) = length(getmin(m))

"""
    getdim(m::AbstractModel)

Returns a 2-tuple of the number of constraints and the number of variables in the model `m`.
"""
getdim(m::AbstractModel) = (getnconstraints(m), getnvars(m))

getmin(m::AbstractModel)= m.box_min
getmin(m::AbstractModel, i) = getmin(m)[i]
isinteger(m::AbstractModel, i) = m.integer[i]

getmax(m::AbstractModel) = m.box_max
getmax(m::AbstractModel, i) = getmax(m)[i]

function getinit(model::AbstractModel)
    _model, _, unflatten = tovecmodel(model)
    return unflatten(getinit(_model))
end

function setmin!(m::AbstractModel, min)
    getmin(m) .= min
    return m
end
function setmin!(m::AbstractModel, i, min)
    getmin(m)[i] = min
    return m
end

function setmax!(m::AbstractModel, max)
    getmax(m) .= max
    return m
end
function setmax!(m::AbstractModel, i, max)
    getmax(m)[i] = max
    return m
end

# Box constraints
function setbox!(m::AbstractModel, minb, maxb)
    getmin(m) .= minb
    getmax(m) .= maxb
    return m
end
function setbox!(m::AbstractModel, i, minb, maxb)
    getmin(m)[i] = minb
    getmax(m)[i] = maxb
    return m
end
function setinteger!(m::AbstractModel, i, integer)
    m.integer[i] = integer
    return m
end

function addvar!(m::Model, lb, ub; init = deepcopy(lb), integer = false)
    push!(getmin(m), lb)
    push!(getmax(m), ub)
    push!(m.init, init)
    push!(m.integer, integer)
    return m
end
function addvar!(m::Model, lb::Vector, ub::Vector; init = deepcopy(lb), integer = falses(length(lb)))
    append!(getmin(m), lb)
    append!(getmax(m), ub)
    append!(m.init, init)
    append!(m.integer, integer)
    return m
end

"""
    set_objective!(m::AbstractModel, f::Function)

Sets the objective of the moodel `m` to the function `f`. `f` must return a scalar.
"""
function set_objective!(m::AbstractModel, f::Function)
    m.objective = Objective(f)
    return m
end

function add_ineq_constraint!(m::AbstractModel, f::Function, s = 0.0; dim = length(flatten(f(getinit(m)))[1]))
    return add_ineq_constraint!(m, FunctionWrapper(f, dim), s)
end
function add_ineq_constraint!(m::AbstractModel, f::AbstractFunction, s = 0.0)
    return add_ineq_constraint!(m, IneqConstraint(f, s))
end
function add_ineq_constraint!(m::AbstractModel, f::IneqConstraint)
    push!(m.ineq_constraints.fs, f)
    return m
end
function add_ineq_constraint!(m::AbstractModel, fs::Vector{<:IneqConstraint})
    append!(m.ineq_constraints.fs, fs)
    return m
end

function add_eq_constraint!(m::AbstractModel, f::Function, s = 0.0; dim = length(flatten(f(getinit(m)))[1]))
    return add_eq_constraint!(m, FunctionWrapper(f, dim), s)
end
function add_eq_constraint!(m::AbstractModel, f::AbstractFunction, s = 0.0)
    return add_eq_constraint!(m, EqConstraint(f, s))
end
function add_eq_constraint!(m::AbstractModel, f::EqConstraint)
    push!(m.eq_constraints.fs, f)
    return m
end
function add_eq_constraint!(m::AbstractModel, fs::Vector{<:EqConstraint})
    append!(m.eq_constraints.fs, fs)
    return m
end
