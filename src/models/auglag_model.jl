# Objective

"""
    AugLagObj

The objective function of the augmented Lagrangian model, `model`. The following are its fields:
 - `model`: the original model optimized
 - `x`: the current primal solution
 - `λ`: the current dual solution
 - `quadweight`: the current quadratic penalty
 - `f`: the current original objective value
 - `g`: the current constraint function value
"""
@params struct AugLagObj <: AbstractFunction
    model::Model
    x::AbstractVector # primal solution
    λ::AbstractVector # dual solution
    quadweight::Base.RefValue # quadratic penalty factor
    f::Base.RefValue
    g::AbstractVector
end

"""
    AugLagObj(model::Model)

Constructs the objective function from the model, `model`.
"""
function AugLagObj(
    model::Model;
    linweights = ones(getnineqconstraints(model)),
    quadweight = 1e-6,
    kwarg...,
)
    x = copy(getmin(model))
    f = Ref(Inf)
    g = similar(linweights)
    return AugLagObj(model, x, linweights, Ref(quadweight), f, g)
end

getdim(::AugLagObj) = 1

"""
    (obj::AugLagObj)(x::AbstractVector, λ::AbstractVector)

Evaluates the augmented Lagrangian at the primal solution `x` and the dual solution `λ`.
"""
function (obj::AugLagObj)(x::AbstractVector, λ::AbstractVector; penalise = true)
    if debugging[]
        #@show x
    end
    parent = getparent(obj)
    origobjval = getobjective(parent)(x)
    origconstrval = getineqconstraints(parent)(x)
    savexλ!(obj, x, λ)
    savefg!(obj, origobjval, origconstrval)

    linepenalty = dot(λ, origconstrval)
    quadpenalty = penalise ? NonNegSumOfSquares(obj.quadweight[])(origconstrval) : zero(linepenalty)

    return origobjval + linepenalty + quadpenalty
end

"""
    (obj::AugLagObj)(primaloptimizer::Function, λ::AbstractVector)

Evaluates the augmented Lagrangian at the dual solution `λ` where `primaloptimizer(λ)` returns a tuple of:
 - The optimal primal solution `x`,
 - The original objective value at the optimal `x` and the current `λ`, and
 - The original constraint values at the optimal `x` and the current `λ`,
"""
function (obj::AugLagObj)(primaloptimizer::Function, λ::AbstractVector; penalise = true)
    parent = getparent(obj)
    x, origobjval, origconstrval = nogradcall(primaloptimizer, λ)
    savexλ!(obj, x, λ)
    savefg!(obj, origobjval, origconstrval)

    linepenalty = dot(λ, origconstrval)
    quadpenalty = penalise ? NonNegSumOfSquares(obj.quadweight[])(origconstrval) : zero(linepenalty)

    return origobjval + linepenalty + quadpenalty
end

nogradcall(f::Function, arg) = f(arg)
ChainRulesCore.@non_differentiable nogradcall(f::Function, arg)

function savefg!(obj::AugLagObj, f::Real, g::AbstractVector)
    setorigobjval!(obj, f)
    setorigconstrval!(obj, g)
    return obj
end
ChainRulesCore.@non_differentiable savefg!(obj::AugLagObj, f::Real, g::AbstractVector)

function savexλ!(obj::AugLagObj, x::AbstractVector, λ::AbstractVector)
    setx!(obj, x)
    setλ!(obj, λ)
    return obj
end
ChainRulesCore.@non_differentiable savexλ!(obj::AugLagObj, x::AbstractVector, λ::AbstractVector)

getparent(f::AugLagObj) = f.model

# To be minimized
#getprimalobjective(f::AugLagObj) = FunctionWrapper(x -> (out = f(x, f.λ); @show out; out), 1)
getprimalobjective(f::AugLagObj; kwargs...) = FunctionWrapper(x -> f(x, f.λ; kwargs...), 1)

# To be maximized
getdualobjective(f::AugLagObj, primaloptimizer::Function) = FunctionWrapper(λ -> begin
    if debugging[]
        #@show λ
    end
    f(primaloptimizer, λ)
end, 1)

getorigobjval(f::AugLagObj) = f.f[]

getorigconstrval(f::AugLagObj) = f.g

getlinweights(f::AugLagObj) = f.λ

getquadweight(f::AugLagObj) = f.quadweight[]

getx(f::AugLagObj) = f.x

getλ(f::AugLagObj) = f.λ

function setlinweights!(f::AugLagObj, λ::AbstractVector)
    getlinweights(f) .= λ
    return f
end

function setquadweight!(f::AugLagObj, quadweight::Number)
    f.quadweight[] = quadweight
    return f
end

function setx!(f::AugLagObj, x::AbstractVector)
    f.x .= x
    return f
end

function setλ!(f::AugLagObj, λ::AbstractVector)
    f.λ .= λ
    return f
end

function setorigobjval!(f::AugLagObj, val::Real)
    f.f[] = val
    return f
end

function setorigconstrval!(f::AugLagObj, g::AbstractVector)
    f.g .= g
    return f
end

# Model

@params struct AugLagModel <: AbstractModel
    parent::Model
    objective::AugLagObj
end
function AugLagModel(
    model::Model;
    kwargs...,
)
    return AugLagModel(model, AugLagObj(model; kwargs...,))
end

getparent(model::AugLagModel) = model.parent

getmin(model::AugLagModel) = getmin(getparent(model))

getmax(model::AugLagModel) = getmax(getparent(model))

getobjective(model::AugLagModel) = model.objective

getineqconstraints(::AugLagModel) = throw("`getineqconstraints` is not defined for `AugLagModel`.")

geteqconstraints(::AugLagModel) = throw("`geteqconstraints` is not defined for `AugLagModel`.")

getobjectiveconstraints(::AugLagModel) = throw("`getobjectiveconstraints` is not defined for `AugLagModel`.")

getlinweights(model::AugLagModel) = getlinweights(getobjective(model))

getx(model::AugLagModel) = getx(getobjective(model))

getλ(model::AugLagModel) = getλ(getobjective(model))

getorigobjval(model::AugLagModel) = getorigobjval(getobjective(model))

getorigconstrval(model::AugLagModel) = getorigconstrval(getobjective(model))

function setlinweights!(model::AugLagModel, λ::AbstractVector)
    setlinweights!(getobjective(model), λ)
    return model
end

getquadweight(model::AugLagModel) = getquadweight(getobjective(model))

function setquadweight!(model::AugLagModel, quadweight::Number)
    setquadweight!(getobjective(model), quadweight)
    return model
end

function setx!(model::AugLagModel, x::AbstractVector)
    setx!(getobjective(model), x)
    return model
end

function setλ!(model::AugLagModel, λ::AbstractVector)
    setλ!(getobjective(model), λ)
    return model
end

function setorigobjval!(model::AugLagModel, val::Real)
    setorigobjval!(getobjective(model), val)
    return model
end

function setorigconstrval!(model::AugLagModel, g::AbstractVector)
    setorigconstrval!(getobjective(model), g)
    return model
end

function getresiduals(solution::Solution, model::AugLagModel, ::GenericCriteria)
    @unpack prevx, x, prevf, f, g = solution
    Δx = maximum(abs(x[j] - prevx[j]) for j in 1:length(x))
    Δf = abs(f - prevf)
    infeas = max(0, maximum(g))
    return Δx, Δf, infeas
end
function getresiduals(solution::Solution, model::AugLagModel, ::KKTCriteria)
    @unpack g, λ, x = solution
    xmin, xmax = getmin(model), getmax(model)
    T = eltype(x)
    auglag = getobjective(model)
    _, ∇L = value_gradient(getprimalobjective(auglag, penalise = true), solution.x)
    res = maximum(1:length(x)) do j
        temp = ∇L[j]
        if xmin[j] >= x[j]
            return abs(min(0, temp))
        elseif x[j] >= xmax[j]
            return max(0, temp)
        else
            return abs(temp)
        end
    end
    if debugging[]
        #@show λ, g
    end
    res = max(res, maximum(abs.(λ .* g)))
    if debugging[]
        @show maximum(abs, g)
        @show maximum(abs, λ)
        @show maximum(x)
    end
    infeas = max(maximum(g), 0)
    return res, infeas
end
function getresiduals(solution::Solution, model::AugLagModel, ::IpoptCriteria)
    @unpack g, λ, x = solution
    if debugging[]
        #@show x
    end
    xmin, xmax = getmin(model), getmax(model)
    T = eltype(x)
    n, m, s = length(x), length(λ), zero(T)
    auglag = getobjective(model)
    _, ∇L = value_gradient(getprimalobjective(auglag, penalise = true), solution.x)
    res = maximum(1:n) do j
        temp = ∇L[j]
        if xmin[j] >= x[j]
            dj = temp
            s += max(dj, 0)
            return abs(min(0, dj))
        elseif x[j] >= xmax[j]
            yj = -temp
            s += max(yj, 0)
            return abs(min(0, yj))
        else
            return abs(temp)
        end
    end
    sd = max(100, (sum(abs, λ) + s) / (n + m)) / 100
    res = max(res, maximum(abs.(λ .* g)))
    res = res / sd
    infeas = max(maximum(g), 0)
    if debugging[]
        println("Agg infeas = ", infeas)
    end
    return res, infeas
end
