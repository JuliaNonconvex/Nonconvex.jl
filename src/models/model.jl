abstract type AbstractModel end

"""
```
struct Model <: AbstractModel
    objective::Objective
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::AbstractVector
    box_max::AbstractVector
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
end

"""
    Model()

Constructs an empty model. The decision variables are assumed to be of type `Vector{Float64}`.
"""
Model() = Model(Float64, nothing)

"""
    Model(f)

Constructs an empty model with objective function `f`. `f` can be an instance of `Base.Function` but must return a number,  or it can be an intance of [`Objective`](@ref). The decision variables are assumed to be of type `Vector{Float64}`.
"""
Model(f::Function) = Model(Float64, f)

"""
    Model(::Type{T}, f::Union{Nothing, Function}) where {T}

Constructs an empty model with objective function `f` and decision variable value type `Vector{T}`. `f` can be an instance of `Base.Function` but must return a number,  or it can be an intance of [`Objective`](@ref). `f` can also be `nothing` in which case, the objective function can be defined later using [`set_objective!`](@ref).
"""
Model(::Type{T}, f::Function) where {T} = Model(T, Objective(f))
function Model(::Type{T}, obj::Union{Nothing, Objective}) where {T}
    return Model(obj, VectorOfFunctions(EqConstraint[]), VectorOfFunctions(IneqConstraint[]), T[], T[])
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

get_objective_multiple(model::Model) = getobjective(model).multiple[]

function set_objective_multiple!(model::Model, m::Real)
    getobjective(model).multiple[] = m
    return model
end

getmin(m::Model)= m.box_min
getmin(m::AbstractModel, i::Integer) = getmin(m)[i]

getmax(m::Model) = m.box_max
getmax(m::AbstractModel, i::Integer) = getmax(m)[i]

function getinit(m::AbstractModel)
    ma = getmax(m)
    mi = getmin(m)
    return map(1:length(mi)) do i
        _ma = ma[i]
        _mi = mi[i]
        _ma == Inf && _mi == -Inf && return 0.0
        _ma == Inf && return _mi + 1.0
        _mi == -Inf && return _ma - 1.0
        return (_ma + _mi) / 2
    end
end

function setmin!(m::AbstractModel, min)
    getmin(m) .= min
    return m
end
function setmin!(m::AbstractModel, i::Integer, min)
    getmin(m)[i] = min
    return m
end

function setmax!(m::AbstractModel, max)
    getmax(m) .= max
    return m
end
function setmax!(m::AbstractModel, i::Integer, max)
    getmax(m)[i] = max
    return m
end

# Box constraints
function setbox!(m::AbstractModel, minb::T, maxb::T) where {T}
    getmin(m) .= minb
    getmax(m) .= maxb
    return m
end
function setbox!(m::AbstractModel, i::Integer, minb, maxb)
    getmin(m)[i] = minb
    getmax(m)[i] = maxb
    return m
end

function addvar!(m::Model, lb, ub)
    append!(getmin(m), lb)
    append!(getmax(m), ub)
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

function add_ineq_constraint!(m::AbstractModel, f::Function, s = 0.0; dim = 1)
    return add_ineq_constraint!(m, FunctionWrapper(f, dim), s)
end
function add_ineq_constraint!(m::AbstractModel, f::AbstractFunction, s = 0.0)
    return add_ineq_constraint!(m, IneqConstraint(f, s))
end
function add_ineq_constraint!(m::Model, f::IneqConstraint)
    push!(m.ineq_constraints.fs, f)
    return m
end
function add_ineq_constraint!(m::Model, fs::Vector{<:IneqConstraint})
    append!(m.ineq_constraints.fs, fs)
    return m
end

function add_eq_constraint!(m::AbstractModel, f::Function, s = 0.0; dim = 1)
    return add_eq_constraint!(m, FunctionWrapper(f, dim), s)
end
function add_eq_constraint!(m::AbstractModel, f::AbstractFunction, s = 0.0)
    return add_eq_constraint!(m, EqConstraint(f, s))
end
function add_eq_constraint!(m::Model, f::EqConstraint)
    push!(m.eq_constraints.fs, f)
    return m
end
function add_eq_constraint!(m::Model, fs::Vector{<:EqConstraint})
    append!(m.eq_constraints.fs, fs)
    return m
end
