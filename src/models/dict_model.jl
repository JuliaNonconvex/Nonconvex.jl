mutable struct DictModel <: AbstractModel
    objective::Union{Nothing, Objective}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    sd_function::Union{AbstractFunction, Nothing}
    box_min::OrderedDict
    box_max::OrderedDict
    init::OrderedDict
    integer::OrderedDict
end

function DictModel(f = nothing)
    return DictModel(
        Objective(f), VectorOfFunctions(EqConstraint[]),
        VectorOfFunctions(IneqConstraint[]),
        nothing,
        OrderedDict(), OrderedDict(),
        OrderedDict(), OrderedDict()
    )
end

function addvar!(m::DictModel, k::Union{Symbol, String}, lb, ub; init = deepcopy(lb), integer = false)
    getmin(m)[k] = lb
    getmax(m)[k] = ub
    m.init[k] = init
    m.integer[k] = integer
    return m
end
