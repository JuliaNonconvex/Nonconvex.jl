"""
```
struct MMADualObj{TApprox <: AbstractMMAApprox} <: AbstractFunction
    model::MMAApproxModel{TApprox}
    x::AbstractVector # primal solution
    approxfg::AbstractVector # approximate [f(x); g(x)]
end
```

The dual objective function of the MMA model, `model`. The approximate objective and constraint values are cached in `approxfg` when the dual objective is called. For every input `λ`, the optimal primal solution is computed and cached in the field `x`.
"""
@params struct MMADualObj{TApprox <: AbstractMMAApprox} <: AbstractFunction
    model::MMAApproxModel{TApprox}
    x::AbstractVector # primal solution
    approxfg::AbstractVector # approximate [f(x); g(x)]
end

"""
    MMADualObj(model::MMAApproxModel)

Constructs the dual objective function from the approximate MMA model, `model`.
"""
function MMADualObj(model::MMAApproxModel{<:MMAApprox})
    xk = getxk(model.approx_objective_ineq_constraints)
    x = copy(xk)
    approxfg = copy(getfk(model.approx_objective_ineq_constraints))
    return MMADualObj(model, x, approxfg)
end
function MMADualObj(model::MMAApproxModel{<:XMMAApprox})
    xk = getxk(model.approx_objective_ineq_constraints)
    yk = model.approx_objective_ineq_constraints.y
    zk = model.approx_objective_ineq_constraints.z[]
    x = [xk; yk; zk]
    approxfg = copy(getfk(model.approx_objective_ineq_constraints))
    return MMADualObj(model, x, approxfg)
end

getdim(::MMADualObj) = 1

# To be maximized
"""
    (f::MMADualObj)(λ)

Evaluates the dual objective value at the dual solution `λ`. The optimal primal solution is computed and used automatically.
"""
function (f::MMADualObj)(λ)
    x = optimizeprimal!(f, λ)
    return f(x, λ)
end

"""
    (f::MMADualObj)(x, λ)

Evaluates the dual objective value at the primal solution `x` and the dual solution `λ`.
"""
function (f::MMADualObj)(x, λ)
    temp = f.model.approx_objective_ineq_constraints(x)
    return temp[1] + dot(temp[2:end], λ)
end

"""
    optimizeprimal!(f::MMADualObj, λ::AbstractVector)

Finds and returns the primal minimizer of the Lagrangian function for a given dual solution `λ`.
"""
function optimizeprimal!(f::MMADualObj{TA}, λ::AbstractVector) where {TA}
    @unpack model, x = f
    @unpack p, q, l, u, r, ρ = getmmaapprox(model.approx_objective_ineq_constraints)
    min = getmin(model)
    max = getmax(model)
    T = eltype(x)

    map!(view(x, 1:size(p, 2)), 1:size(p, 2)) do j
        @views t1 = p[1, j] + λ' * p[2:end, j]
        @views t2 = q[1, j] + λ' * q[2:end, j]
        temp1 = t1 / (u[j] - min[j])^2 - t2 / (min[j] - l[j])^2
        temp2 = t1 / (u[j] - max[j])^2 - t2 / (max[j] - l[j])^2
        if isfinite(min[j]) && temp1 >= 0
            return min[j]
        elseif isfinite(max[j]) && temp2 <= 0
            return max[j]
        else
            return (sqrt(t1) * l[j] + sqrt(t2) * u[j]) / (sqrt(t1) + sqrt(t2))
        end
    end

    if TA <: XMMAApprox
        @unpack a, c, d = model.approx_objective_ineq_constraints
        map!(view(x, size(p, 2) + 1 : size(p, 2) + length(c)), 1:length(c)) do i
            temp = c[i] - λ[i]
            if temp >= 0
                return zero(T)
            else
                return -temp / d[i]
            end
        end
        temp = a[1] - dot(a[2:end], λ)
        ε = 1e-4
        if temp >= 0
            x[end] = zero(T)
        else
            x[end] = -temp / (2 * ε)
        end
    end

    return x
end
ChainRulesCore.@non_differentiable optimizeprimal!(f::MMADualObj, λ::AbstractVector)

"""
    getoptimalx(f::MMADualObj)
    
Returns the cached optimal primal solution.
"""
getoptimalx(f::MMADualObj) = f.x[1:getnvars(getparent(f.model))]

"""
    getapproxfg(dualobj::MMADualObj)

Returns the cached approximate objective and constraint values.
"""
getapproxfg(dualobj::MMADualObj) = getapproxfg(dualobj.model)

"""
```
struct MMADualModel <: AbstractModel
    parent::MMAApproxModel
    obj::MMADualObj
end
```

The dual model with a primal model `parent` and dual objective `obj`.
"""
@params struct MMADualModel <: AbstractModel
    parent::MMAApproxModel
    obj::MMADualObj
end

"""
    MMADualModel(model::MMAApproxModel)

Constructs an instance of `MMADualModel` for the approximate MMA model, `model`.
"""
function MMADualModel(model::MMAApproxModel)
    obj = MMADualObj(model)
    return MMADualModel(model, obj)
end

"""
    getoptimalx(m::MMADualModel)

Returns the last cached optimal primal solution.
"""
getoptimalx(m::MMADualModel) = getoptimalx(m.obj)

getparent(m::MMADualModel) = m.parent

@doc docσ
getσ(m::MMADualModel) = getσ(getparent(m))

@doc docσ
setσ!(m::MMADualModel, σ) = setσ!(getparent(m), σ)

@doc docρ
getρ(m::MMADualModel) = getρ(getparent(m))

@doc docρ
setρ!(m::MMADualModel, ρ) = setρ!(getparent(m), ρ)

@doc docupdateapprox!
function updateapprox!(m::MMADualModel, args...)
    updateapprox!(getparent(m), args...)
    return m
end

getobjective(m::MMADualModel) = m.obj

@doc docxk
getxk(m::MMADualModel) = getxk(getparent(m))

@doc docfk
getfk(m::MMADualModel) = getfk(getparent(m))

@doc doc∇fk
get∇fk(m::MMADualModel) = get∇fk(getparent(m))
