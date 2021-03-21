struct PercivalAlg end

@params struct PercivalOptions
    nt::NamedTuple
end
function PercivalOptions(; first_order = true, memory = 5, kwargs...)
    return PercivalOptions(
        (;first_order = first_order, memory = memory, kwargs...),
    )
end

@params mutable struct PercivalWorkspace <: Workspace
    model::Model
    problem::Percival.NLPModels.AbstractNLPModel
    x0::AbstractVector
    options::PercivalOptions
end
function PercivalWorkspace(
    model::Model, x0::AbstractVector = getinit(model);
    options = PercivalOptions(), kwargs...,
)
    problem = getpercival_problem(model, x0)
    return PercivalWorkspace(model, problem, x0, options)
end
@params struct PercivalResult
    minimizer
    minimum
    problem
    result
end

function optimize!(workspace::PercivalWorkspace)
    @unpack problem, options, x0 = workspace
    foreach(keys(options.nt)) do k
        if k != :first_order && k != :memory
            v = options.nt[k]
            setproperty!(problem, k, v)
        end
    end
    if options.nt.first_order
        qnlp = Percival.LBFGSModel(
            Percival.SlackModel(problem);
            mem = options.nt.memory,
        )
        result = percival(qnlp)
        result.solution = result.solution[1:length(x0)]
    else
        result = percival(problem)
    end
    return PercivalResult(
        copy(result.solution), result.objective, problem, result,
    )
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
    return getpercival_problem(
        getobjective(model),
        ineq,
        eq,
        x0,
        getmin(model),
        getmax(model),
    )
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
