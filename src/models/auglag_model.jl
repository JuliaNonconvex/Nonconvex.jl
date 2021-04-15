# Objective

"""
    AugLag2Obj

The objective function of the augmented Lagrangian model, `model`. The following are its fields:
 - `model`: the original model optimized
 - `x`: the current primal solution
 - `λ`: the current dual solution
 - `quadweight`: the current quadratic penalty
 - `f`: the current original objective value
 - `g`: the current constraint function value
"""
@params struct AugLag2Obj <: AbstractFunction
    model::Model
    x::AbstractVector # primal solution
    λ::AbstractVector # dual solution
    quadweight::Base.RefValue # quadratic penalty factor
    f::Base.RefValue
    g::AbstractVector
end

"""
    AugLag2Obj(model::Model)

Constructs the objective function from the model, `model`.
"""
function AugLag2Obj(
    model::Model;
    linweights = ones(getnineqconstraints(model)),
    quadweight = 1e-6,
    kwarg...,
)
    x = copy(getmin(model))
    f = Ref(Inf)
    g = similar(linweights)
    return AugLag2Obj(model, x, linweights, Ref(quadweight), f, g)
end

getdim(::AugLag2Obj) = 1

"""
    (obj::AugLag2Obj)(x::AbstractVector, λ::AbstractVector)

Evaluates the augmented Lagrangian at the primal solution `x` and the dual solution `λ`.
"""
function (obj::AugLag2Obj)(x::AbstractVector, λ::AbstractVector; penalise = true)
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
    (obj::AugLag2Obj)(primaloptimizer::Function, λ::AbstractVector)

Evaluates the augmented Lagrangian at the dual solution `λ` where `primaloptimizer(λ)` returns a tuple of:
 - The optimal primal solution `x`,
 - The original objective value at the optimal `x` and the current `λ`, and
 - The original constraint values at the optimal `x` and the current `λ`,
"""
function (obj::AugLag2Obj)(primaloptimizer::Function, λ::AbstractVector; penalise = true)
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

function savefg!(obj::AugLag2Obj, f::Real, g::AbstractVector)
    setorigobjval!(obj, f)
    setorigconstrval!(obj, g)
    return obj
end
ChainRulesCore.@non_differentiable savefg!(obj::AugLag2Obj, f::Real, g::AbstractVector)

function savexλ!(obj::AugLag2Obj, x::AbstractVector, λ::AbstractVector)
    setx!(obj, x)
    setλ!(obj, λ)
    return obj
end
ChainRulesCore.@non_differentiable savexλ!(obj::AugLag2Obj, x::AbstractVector, λ::AbstractVector)

getparent(f::AugLag2Obj) = f.model

# To be minimized
#getprimalobjective(f::AugLag2Obj) = FunctionWrapper(x -> (out = f(x, f.λ); @show out; out), 1)
getprimalobjective(f::AugLag2Obj; kwargs...) = FunctionWrapper(x -> f(x, f.λ; kwargs...), 1)

# To be maximized
getdualobjective(f::AugLag2Obj, primaloptimizer::Function) = FunctionWrapper(λ -> begin
    if debugging[]
        #@show λ
    end
    f(primaloptimizer, λ)
end, 1)

getorigobjval(f::AugLag2Obj) = f.f[]

getorigconstrval(f::AugLag2Obj) = f.g

getlinweights(f::AugLag2Obj) = f.λ

getquadweight(f::AugLag2Obj) = f.quadweight[]

getx(f::AugLag2Obj) = f.x

getλ(f::AugLag2Obj) = f.λ

function setlinweights!(f::AugLag2Obj, λ::AbstractVector)
    getlinweights(f) .= λ
    return f
end

function setquadweight!(f::AugLag2Obj, quadweight::Number)
    f.quadweight[] = quadweight
    return f
end

function setx!(f::AugLag2Obj, x::AbstractVector)
    f.x .= x
    return f
end

function setλ!(f::AugLag2Obj, λ::AbstractVector)
    f.λ .= λ
    return f
end

function setorigobjval!(f::AugLag2Obj, val::Real)
    f.f[] = val
    return f
end

function setorigconstrval!(f::AugLag2Obj, g::AbstractVector)
    f.g .= g
    return f
end

# Model

@params struct AugLag2Model <: AbstractModel
    parent::Model
    objective::AugLag2Obj
end
function AugLag2Model(
    model::Model;
    kwargs...,
)
    return AugLag2Model(model, AugLag2Obj(model; kwargs...,))
end

getparent(model::AugLag2Model) = model.parent

getmin(model::AugLag2Model) = getmin(getparent(model))

getmax(model::AugLag2Model) = getmax(getparent(model))

getobjective(model::AugLag2Model) = model.objective

getineqconstraints(::AugLag2Model) = throw("`getineqconstraints` is not defined for `AugLag2Model`.")

geteqconstraints(::AugLag2Model) = throw("`geteqconstraints` is not defined for `AugLag2Model`.")

getobjectiveconstraints(::AugLag2Model) = throw("`getobjectiveconstraints` is not defined for `AugLag2Model`.")

getlinweights(model::AugLag2Model) = getlinweights(getobjective(model))

getx(model::AugLag2Model) = getx(getobjective(model))

getλ(model::AugLag2Model) = getλ(getobjective(model))

getorigobjval(model::AugLag2Model) = getorigobjval(getobjective(model))

getorigconstrval(model::AugLag2Model) = getorigconstrval(getobjective(model))

function setlinweights!(model::AugLag2Model, λ::AbstractVector)
    setlinweights!(getobjective(model), λ)
    return model
end

getquadweight(model::AugLag2Model) = getquadweight(getobjective(model))

function setquadweight!(model::AugLag2Model, quadweight::Number)
    setquadweight!(getobjective(model), quadweight)
    return model
end

function setx!(model::AugLag2Model, x::AbstractVector)
    setx!(getobjective(model), x)
    return model
end

function setλ!(model::AugLag2Model, λ::AbstractVector)
    setλ!(getobjective(model), λ)
    return model
end

function setorigobjval!(model::AugLag2Model, val::Real)
    setorigobjval!(getobjective(model), val)
    return model
end

function setorigconstrval!(model::AugLag2Model, g::AbstractVector)
    setorigconstrval!(getobjective(model), g)
    return model
end

function getresiduals(solution::Solution, model::AugLag2Model, ::GenericCriteria)
    @unpack prevx, x, prevf, f, g = solution
    Δx = maximum(abs(x[j] - prevx[j]) for j in 1:length(x))
    Δf = abs(f - prevf)
    infeas = max(0, maximum(g))
    return Δx, Δf, infeas
end
function getresiduals(solution::Solution, model::AugLag2Model, ::KKTCriteria)
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
function getresiduals(solution::Solution, model::AugLag2Model, ::IpoptCriteria)
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
