"""
```
struct MMALagModel <: AbstractModel
    parent::VecModel
    quadweight::Base.RefValue{<:Real}
    linweights::AbstractVector
    origobjval::Base.RefValue{<:Real}
    origconstrval::AbstractVector
end
```

A model with aggregated block constraints. Any block constraint in `parent` is aggregated using the corresponding weights in `linweights`. Non-block constraints are treated as block constraints of dimension 1. Let `violations` be the sum of constraint violations at a point `x`. The objective value of an instance of `MMALagModel` is the original objective value of its `parent` plus `quadweight[] * sum(max.(violations, 0)^2)`. The Lagrangian function of the new problem is the augmented Lagrangian function of the original problem. `origobjval` and `origconstrval` store the original objective and constraint values respectively when evaluating the aggregated objective and constraints.
"""
@params struct MMALagModel <: AbstractModel
    parent::VecModel
    quadweight::Base.RefValue{<:Real}
    linweights::AbstractVector
    origobjval::Base.RefValue{<:Real}
    origconstrval::AbstractVector
end
function MMALagModel(
    model::VecModel;
    quadweight = 1e-5,
    linweights = ones(2*getdim(getineqconstraints(model))),
    kwargs...,
)
    return MMALagModel(
        model,
        Ref(quadweight),
        linweights,
        Ref(zero(quadweight)),
        zeros(getdim(getineqconstraints(model))),
    )
end
MMALagModel(args...; kwargs...) = MMALagModel(Model(args...); kwargs...)

getparent(m::MMALagModel) = m.parent

function getobjectiveconstraints(m::MMALagModel)
    parent = getparent(m)
    quadpenalty = NonNegSumOfSquares(m.quadweight[])

    dims = getdim.(getineqconstraints(getparent(m)).fs)
    M = sum(dims)
    nblocks = length(dims)
    func = x -> begin
        origobjval = getobjective(parent)(x)
        origconstrval = getineqconstraints(parent)(x)
        setvals!(m, origobjval / get_objective_multiple(parent), origconstrval)
        
        obj = origobjval + quadpenalty(origconstrval)
        start = 1
        constr1 = map(1:nblocks) do i
            inds = start:(start+dims[i]-1)
            start += dims[i]
            return @views WeightedSum(m.linweights[inds])(origconstrval[inds])
        end
        constr2 = map(1:nblocks) do i
            inds = start:(start+dims[i]-1)
            start += dims[i]
            return @views WeightedSum(m.linweights[inds])(origconstrval[inds .- M])
        end
        @assert start == 2*length(origconstrval) + 1
        return [obj; constr1; constr2]
    end
    return FunctionWrapper(func, 2nblocks + 1)
end

function getobjective(m::MMALagModel)
    objf(x) = getobjectiveconstraints(m)(x)[1]
    return Objective(objf)
end

function getineqconstraint(m::MMALagModel, i::Integer)
    throw("`getineqconstraint` is not defined for `MMALagModel`.")
end

function geteqconstraint(m::MMALagModel, i::Integer)
    throw("`geteqconstraint` is not defined for `MMALagModel`.")
end

function getineqconstraints(m::MMALagModel)
    dim = getdim(getobjectiveconstraints(m)) - 1
    constrf(x) = getobjectiveconstraints(m)(x)[2:end]
    return FunctionWrapper(constrf, dim)
end

function geteqconstraints(m::MMALagModel)
    throw("`geteqconstraint` is not defined for `MMALagModel`.")
end

get_objective_multiple(model::MMALagModel) = get_objective_multiple(getparent(model))

function set_objective_multiple!(model::MMALagModel, m::Real)
    set_objective_multiple!(getparent(model), m)
    return model
end

getmin(m::MMALagModel)= getmin(m.parent)

getmax(m::MMALagModel) = getmax(m.parent)

function addvar!(m::MMALagModel, lb, ub)
    addvar!(getparent(m), lb, ub)
    return m
end

function add_ineq_constraint!(m::MMALagModel, f::IneqConstraint)
    add_ineq_constraint!(getparent(m), f)
    append!(m.linweights, ones(2*getdim(f)))
    append!(m.origconstrval, zeros(getdim(f)))
    return m
end

function add_ineq_constraint!(m::MMALagModel, fs::Vector{<:IneqConstraint})
    add_ineq_constraint!(getparent(m), fs)
    append!(m.linweights, ones(2*getdim(VectorOfFunctions(fs))))
    append!(m.origconstrval, zeros(getdim(VectorOfFunctions(fs))))
    return m
end

"""
    getquadweight(m::MMALagModel)

Returns the quadratic penalty of `m`.
"""
getquadweight(m::MMALagModel) = m.quadweight[]

"""
    getlinweights(m::MMALagModel)

Returns the linear weights of `m`. These are the weights used to aggregate the block constraints in the `MMALag` algorithm.
"""
getlinweights(m::MMALagModel) = m.linweights

"""
    getorigobjval(m::MMALagModel)

Returns the last cached objective value.
"""
getorigobjval(m::MMALagModel) = m.origobjval[]

"""
    getorigconstrval(m::MMALagModel)

Returns the last cached constraint values.
"""
getorigconstrval(m::MMALagModel) = m.origconstrval

"""
    setvals!(m::MMALagModel, obj, constr)

Caches the objective and constraint values `obj` and `constr` in `m`.
"""
function setvals!(m::MMALagModel, obj, constr)
    m.origobjval[] = obj
    m.origconstrval .= constr
    return m
end
ChainRulesCore.@non_differentiable setvals!(m::MMALagModel, obj, constr)

"""
    setquadweight!(m::MMALagModel, w)

Sets the quadratic penalty in `m` to `w`.
"""
function setquadweight!(m::MMALagModel, w)
    m.quadweight[] = w
    return m
end

"""
    setlinweights!(m::MMALagModel, λ)

Sets linear weights in `m` to `λ`.
"""
function setlinweights!(m::MMALagModel, λ)
    m.linweights .= λ
    return m
end
