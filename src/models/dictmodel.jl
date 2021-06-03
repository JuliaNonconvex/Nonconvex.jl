mutable struct DictModel <: AbstractModel
    objective::Union{Nothing, Function}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::OrderedDict
    box_max::OrderedDict
end
function tomodel(m::DictModel, x0 = getmin(m))
    v, unflatten = flatten(x0)
    return Model(
        Objective(x -> m.objective(unflatten(x))),
        length(m.eq_constraints.fs) != 0 ? VectorOfFunctions(map(m.eq_constraints.fs) do c
            EqConstraint(x -> c.f(unflatten(x)), c.rhs, c.dim)
        end) : VectorOfFunctions(EqConstraint[]),
        length(m.ineq_constraints.fs) != 0 ? VectorOfFunctions(map(m.ineq_constraints.fs) do c
            IneqConstraint(x -> c.f(unflatten(x)), c.rhs, c.dim)
        end) : VectorOfFunctions(IneqConstraint[]),
        flatten(m.box_min)[1],
        flatten(m.box_max)[1],
    ), v, unflatten
end

function DictModel(f = nothing)
    return DictModel(
        f, VectorOfFunctions(EqConstraint[]),
        VectorOfFunctions(IneqConstraint[]),
        OrderedDict(), OrderedDict(),
    )
end

function addvar!(m::DictModel, k::Union{Symbol, String}, lb, ub)
    getmin(m)[k] = lb
    getmax(m)[k] = ub
    return m
end

function optimize(model::DictModel, optimizer::AbstractOptimizer, x0, args...; kwargs...)
    _model, _x0, unflatten = tomodel(model, x0)
    r = optimize(_model, optimizer, _x0, args...; kwargs...)
    return @set r.minimizer = unflatten(r.minimizer)
end

function getinit(model::DictModel)
    _model, _, unflatten = tomodel(model)
    return unflatten(getinit(_model))
end

struct Unflatten{F} <: Function
    unflatten::F
end
(f::Unflatten)(x) = f.unflatten(x)

function _merge(d1, d2)
    _d = OrderedDict(k => zero(v) for (k, v) in d1)
    return merge(_d, d2)
end

function ChainRulesCore.rrule(un::Unflatten, v)
    x = un.unflatten(v)
    return x, Δ -> (NO_FIELDS, flatten(_merge(x, Δ))[1])
end

flatten(x) = ParameterHandling.flatten(x)
function flatten(d::OrderedDict)
    v, un = flatten(collect(values(d)))
    return v, Unflatten(x -> OrderedDict(collect(keys(d)) .=> un(x)))
end

function ChainRulesCore.rrule(::typeof(flatten), x)
    v, un = flatten(x)
    return (v, un), Δ -> (NO_FIELDS, un(Δ[1]))
end
