mutable struct VecModel{Tv <: AbstractVector} <: AbstractModel
    objective::Union{Nothing, Objective}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    sd_function::Union{MatrixFunctionWrapper, Nothing}
    box_min::Tv
    box_max::Tv
    init::Tv
    integer::BitVector
end
function VecModel(
        objective::Union{Nothing, Objective}, 
        eq_constraints::VectorOfFunctions, 
        ineq_constraints::VectorOfFunctions, 
        box_min, 
        box_max, 
        init, 
        integer::BitVector)
    return VecModel(objective, eq_constraints, ineq_constraints, nothing, box_min, box_max, init, integer)
end

function isfeasible(model::VecModel, x::AbstractVector; ctol = 1e-4)
    return all(getmin(model) .<= x .<= getmax(model)) &&
        all(getineqconstraints(model)(x) .<= ctol) &&
        all(-ctol .<= geteqconstraints(model)(x) .<= ctol)
end

function addvar!(m::VecModel, lb::Real, ub::Real; init::Real = lb, integer = false)
    push!(getmin(m), lb)
    push!(getmax(m), ub)
    push!(m.init, init)
    push!(m.integer, integer)
    return m
end
function addvar!(m::VecModel, lb::Vector{<:Real}, ub::Vector{<:Real}; init::Vector{<:Real} = copy(lb), integer = falses(length(lb)))
    append!(getmin(m), lb)
    append!(getmax(m), ub)
    append!(m.init, init)
    append!(m.integer, integer)
    return m
end

function getinit(m::VecModel)
    ma = getmax(m)
    mi = getmin(m)
    init = m.init
    return map(1:length(mi)) do i
        if isfinite(init[i])
            return init[i]
        else
            _ma = ma[i]
            _mi = mi[i]
            _ma == Inf && _mi == -Inf && return 0.0
            _ma == Inf && return _mi + 1.0
            _mi == -Inf && return _ma - 1.0
            return (_ma + _mi) / 2
        end
    end
end

get_objective_multiple(model::VecModel) = getobjective(model).multiple[]

function set_objective_multiple!(model::VecModel, m)
    getobjective(model).multiple[] = m
    return model
end

function optimize(model::VecModel, optimizer::AbstractOptimizer, x0, args...; kwargs...)
    workspace = Workspace(model, optimizer, copy(x0), args...; kwargs...)
    return optimize!(workspace)
end

function tovecmodel(m::AbstractModel, x0 = getmin(m))
    v, _unflatten = flatten(x0)
    unflatten = Unflatten(x0, _unflatten)
    # obj = Objective(x -> m.objective(unflatten(x)), flags = m.objective.flags)
    # wrapper = MatrixFunctionWrapper(x -> m.sd_function(unflatten(x)), m.sd_function.mat_dim)
    return VecModel(
        # objective
        Objective(x -> m.objective(unflatten(x)), flags = m.objective.flags),
        # eq_constraints
        length(m.eq_constraints.fs) != 0 ? VectorOfFunctions(map(m.eq_constraints.fs) do c
            EqConstraint(x -> maybeflatten(c.f(unflatten(x)))[1], maybeflatten(c.rhs)[1], c.dim, c.flags)
        end) : VectorOfFunctions(EqConstraint[]),
        # ineq_constraints
        length(m.ineq_constraints.fs) != 0 ? VectorOfFunctions(map(m.ineq_constraints.fs) do c
            IneqConstraint(x -> maybeflatten(c.f(unflatten(x)))[1], maybeflatten(c.rhs)[1], c.dim, c.flags)
        end) : VectorOfFunctions(IneqConstraint[]),
        # sd_function
        m.sd_function isa Nothing ? nothing : MatrixFunctionWrapper(x -> m.sd_function(unflatten(x)), m.sd_function.mat_dim), 
        # box_min
        float.(flatten(m.box_min)[1]),
        # box_max
        float.(flatten(m.box_max)[1]),
        # init
        float.(flatten(m.init)[1]),
        # integer
        convert(BitVector, flatten(m.integer)[1]),
    ), float.(v), unflatten
end
