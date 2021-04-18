struct PercivalAlg end
const AugLag = PercivalAlg

@params struct PercivalOptions
    nt::NamedTuple
end
function PercivalOptions(; first_order = true, memory = 5, kwargs...)
    return PercivalOptions(
        (;first_order = first_order,
        memory = memory, inity = ones,
        kwargs...),
    )
end
const AugLagOptions = PercivalOptions

@params mutable struct PercivalWorkspace <: Workspace
    model::Model
    problem::Percival.NLPModels.AbstractNLPModel
    x0::AbstractVector
    options::PercivalOptions
    counter::Base.RefValue{Int}
end
function PercivalWorkspace(
    model::Model, x0::AbstractVector = getinit(model);
    options = PercivalOptions(), kwargs...,
)
    problem, counter = getpercival_problem(model, x0)
    return PercivalWorkspace(model, problem, x0, options, counter)
end
@params struct PercivalResult<:AbstractResult
    minimizer
    minimum
    problem
    result
    fcalls
end

function optimize!(workspace::PercivalWorkspace)
    @unpack problem, options, x0, counter = workspace
    counter[] = 0
    m = getnconstraints(workspace.model)
    result = _percival(problem; options.nt..., inity = options.nt.inity(m))
    result.solution = result.solution[1:length(x0)]
    return PercivalResult(
        copy(result.solution), result.objective, problem, result, counter[],
    )
end

function _percival(nlp;
    μ::Real = eltype(nlp.meta.x0)(10.0),
    max_iter::Int = 1000, max_time::Real = Inf,
    max_eval::Int = 100000, atol::Real = 1e-6,
    rtol::Real = 1e-6, ctol::Real = 1e-6, first_order = true, memory = 5,
    subsolver_logger::Percival.AbstractLogger = Percival.NullLogger(),
    inity = nothing, max_cgiter = 100, subsolver_max_eval = 200, kwargs...,
)
    modifier = m -> NLPModelsModifiers.LBFGSModel(m, mem = memory)
    _kwargs = (
        max_iter = max_iter, max_time = max_time,
        max_eval = max_eval, atol = atol, rtol = rtol,
        subsolver_logger = subsolver_logger,
        subproblem_modifier = first_order ? modifier : identity,
        subsolver_max_eval = subsolver_max_eval,
        subsolver_kwargs = Dict(:max_cgiter => max_cgiter),
    )
    if Percival.unconstrained(nlp) || Percival.bound_constrained(nlp)
        return Percival.percival(
            Val(:tron), nlp; _kwargs...,
        )
    elseif Percival.equality_constrained(nlp)
        return Percival.percival(
            Val(:equ), nlp; inity = inity, _kwargs...,
        )
    else # has inequalities
        return Percival.percival(
            Val(:ineq), nlp; inity = inity, _kwargs...,
        )
    end
end

function Workspace(model::AbstractModel, optimizer::PercivalAlg, args...; kwargs...,)
    return PercivalWorkspace(model, args...; kwargs...)
end

function getpercival_problem(model::Model, x0::AbstractVector)
    eq = if length(model.eq_constraints.fs) == 0
        nothing
    else
        model.eq_constraints
    end
    ineq = if length(model.ineq_constraints.fs) == 0
        nothing
    else
        model.ineq_constraints
    end
    obj = CountingFunction(getobjective(model))
    return getpercival_problem(
        obj,
        ineq,
        eq,
        x0,
        getmin(model),
        getmax(model),
    ), obj.counter
end
function getpercival_problem(obj, ineq_constr, eq_constr, x0, xlb, xub)
    nvars = length(x0)
    if ineq_constr !== nothing
        ineqval = ineq_constr(x0)
        ineq_nconstr = length(ineqval)
    else
        ineqval = Float64[]
        ineq_nconstr = 0
    end
    if eq_constr !== nothing
        eqval = eq_constr(x0)
        eq_nconstr = length(eqval)
    else
        eqval = Float64[]
        eq_nconstr = 0
    end
    c = x -> begin
        if ineq_constr !== nothing
            v1 = ineq_constr(x)
        else
            v1 = eltype(x)[]
        end
        if eq_constr !== nothing
            v2 = eq_constr(x)
        else
            v2 = eltype(x)[]
        end
        return [v1; v2]
    end
    lcon = [fill(-Inf, ineq_nconstr); zeros(eq_nconstr)]
    ucon = zeros(ineq_nconstr + eq_nconstr)
    nlp = ADNLPModels.ADNLPModel(obj, x0, xlb, xub, c, lcon, ucon, adbackend = ADNLPModels.ZygoteAD())
    return nlp
end
