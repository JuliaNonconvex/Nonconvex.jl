"""
```
struct MMAApproxModel{TApprox <: AbstractMMAApprox} <: AbstractModel
    parent::VecModel
    objective_ineq_constraints::AbstractFunction
    approx_objective_ineq_constraints::TApprox
    box_min::AbstractVector
    box_max::AbstractVector
end
```

An approximate restricted model that uses the `MMAApprox` or `XMMAApprox` approximation of the objective and constraint functions. The lower and upper bounds of the approximate model are more restricted than the original model. See [`updatebounds!`](@ref) for more details on the bounds.
- `parent`: the original model.
- `objective_ineq_constraints`: a function that returns the exact objective value and constraint violations as a single vector.
- `approx_objective_ineq_constraints`: a function that returns the approximate objective value and constraint violations as a single vector.
- `box_min`: the restricted lower bounds on the decision variables.
- `box_max`: the restricted upper bounds on the decision variables.
"""
@params struct MMAApproxModel{TApprox <: AbstractMMAApprox} <: AbstractModel
    parent::VecModel
    objective_ineq_constraints::AbstractFunction
    approx_objective_ineq_constraints::TApprox
    box_min::AbstractVector
    box_max::AbstractVector
end

"""
    MMAApproxModel(parent::VecModel, x::AbstractVector; extended = false, kwargs...)

Constructs an MMA approximation of the model `parent` around the point `x`. If `extended` is true, an [`XMMAApprox`](@ref) approximation is used. Otherwise, the default [`MMAApprox`](@ref) approximation is used. `kwargs` can be used to pass additional options to `XMMAApprox`, e.g. setting the coefficients.
"""
function MMAApproxModel(
    parent::VecModel,
    x::AbstractVector;
    extended = false,
    kwargs...,
)
    T = eltype(x)
    σ = map(1:length(x)) do j
        diff = getmax(parent, j) - getmin(parent, j)
        if !isfinite(diff)
            diff = 1000 * one(T)
        end
        return 0.5 * diff
    end
    ρ = zeros(length(getineqconstraints(parent)) + 1)
    obj_constr = getobjectiveconstraints(parent)
    approx_obj_constr = MMAApprox(
        obj_constr,
        x;
        σ = σ,
        ρ = ρ,
    )
    model = MMAApproxModel(
        parent,
        obj_constr,
        extended ? XMMAApprox(approx_obj_constr; kwargs...) : approx_obj_constr,
        copy(getmin(parent)),
        copy(getmax(parent)),
    )
    updateapprox!(model)
    return model
end

getmin(m::MMAApproxModel)= m.box_min
getmax(m::MMAApproxModel) = m.box_max

@doc docρ
getρ(model::MMAApproxModel) = getρ(model.approx_objective_ineq_constraints)

@doc docρ
function setρ!(model::MMAApproxModel, ρ)
    setρ!(model.approx_objective_ineq_constraints, ρ)
    return model
end

@doc docσ
getσ(model::MMAApproxModel) = getσ(model.approx_objective_ineq_constraints)

@doc docσ
function setσ!(model::MMAApproxModel, σ)
    setσ!(model.approx_objective_ineq_constraints, σ)
    return model
end

@doc docxk
getxk(model::MMAApproxModel) = getxk(model.approx_objective_ineq_constraints)

@doc docxk
function setxk!(model::MMAApproxModel, x::AbstractVector)
    setxk!(model.approx_objective_ineq_constraints, x)
    return model
end

@doc docfk
getfk(model::MMAApproxModel) = getfk(model.approx_objective_ineq_constraints)

@doc docfk
function setfk!(model::MMAApproxModel, f)
    setfk!(model.approx_objective_ineq_constraints, f)
    return model
end

@doc doc∇fk
get∇fk(model::MMAApproxModel) = get∇fk(model.approx_objective_ineq_constraints)

@doc doc∇fk
function set∇fk!(model::MMAApproxModel, ∇f::AbstractVecOrMat)
    set∇fk!(model.approx_objective_ineq_constraints, ∇f)
    return model
end

getparent(model::MMAApproxModel) = model.parent

"""
    getapproxfg(model::MMAApproxModel)

Returns the last cached value of the approximate objective and constraint values.
"""
getapproxfg(model::MMAApproxModel) = getapproxfg(model.approx_objective_ineq_constraints)

"""
    updatebounds!(model::MMAApproxModel, x)

Updates the bounds on the decision variables in `model` using the asymptotes of the approximation and the bounds on the original variables. The update formula is described in the original [1987 MMA paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207) and the [2002 paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822).
"""
function updatebounds!(model::MMAApproxModel, x)
    L, U = getasymptotes(model.approx_objective_ineq_constraints)
    getmin(model) .= max.(getmin(getparent(model)), L .+ 0.1 .* (x .- L))
    getmax(model) .= min.(getmax(getparent(model)), U .- 0.1 .* (U .- x))
    return model
end

# Must be called every time σ or ρ are updated
# Doesn't do function evaluations
"""
    updateapprox!(approx::MMAApproxModel)

Updates the MMA approximation in `model`. This is called when `σ`, `ρ` or `xk` change.
"""
function updateapprox!(model::MMAApproxModel)
    updateapprox!(model.approx_objective_ineq_constraints)
    updatebounds!(model, getxk(model))
    return model
end

# Performs a function evaluation at x
"""
    updateapprox!(model::MMAApproxModel, x)

Updates the MMA approximation in `model` making it around the point `x`.
"""
function updateapprox!(model::MMAApproxModel, x::AbstractVector)
    updateapprox!(model.approx_objective_ineq_constraints, x)
    updatebounds!(model, x)
    return model
end

# Doesn't do function evaluations
"""
    updateapprox!(model::MMAApproxModel, x, f, ∇f)

Updates the MMA approximation in `model` making it around the point `x` where `f` is the value of the approximated function at `x` and `∇f` is its gradient if `f` is a number and the Jacobian if `f` is a vector.
"""
function updateapprox!(model::MMAApproxModel, x::AbstractVector, f, ∇f)
    updateapprox!(model.approx_objective_ineq_constraints, x, f, ∇f)
    updatebounds!(model, x)
    return model
end
