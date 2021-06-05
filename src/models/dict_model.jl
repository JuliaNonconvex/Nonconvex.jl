mutable struct DictModel <: AbstractModel
    objective::Union{Nothing, Function}
    eq_constraints::VectorOfFunctions
    ineq_constraints::VectorOfFunctions
    box_min::OrderedDict
    box_max::OrderedDict
    init::OrderedDict
    integer::OrderedDict
end

function DictModel(f = nothing)
    return DictModel(
        f, VectorOfFunctions(EqConstraint[]),
        VectorOfFunctions(IneqConstraint[]),
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

function _getinteger(m::DictModel)
    return convert(BitVector, reduce(vcat, fill.(values(m.integer), length.(getindex.(flatten.(values(m.box_min)), 1)))))
end
