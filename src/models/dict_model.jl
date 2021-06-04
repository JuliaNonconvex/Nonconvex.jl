mutable struct DictModel <: AbstractModel
    objective::Union{Nothing, Function}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::OrderedDict
    box_max::OrderedDict
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
