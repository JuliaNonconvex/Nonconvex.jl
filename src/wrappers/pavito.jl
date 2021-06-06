import .Pavito

@params struct PavitoIpoptCbcOptions
    nt::NamedTuple
    subsolver_options
    first_order::Bool
end
function PavitoIpoptCbcOptions(;
    first_order = true, subsolver_options = IpoptOptions(first_order = first_order), kwargs...,
)
    first_order = !hasproperty(subsolver_options.nt, :hessian_approximation) ||
        subsolver_options.nt.hessian_approximation == "limited-memory"
    return PavitoIpoptCbcOptions((; kwargs...), subsolver_options, first_order)
end

@params mutable struct PavitoIpoptCbcWorkspace <: Workspace
    model::VecModel
    problem::JuMPProblem
    x0::AbstractVector
    integers::AbstractVector{<:Bool}
    options::PavitoIpoptCbcOptions
    counter::Base.RefValue{Int}
end
function PavitoIpoptCbcWorkspace(
    model::VecModel, x0::AbstractVector = getinit(model);
    options = PavitoIpoptCbcOptions(), kwargs...,
)
    integers = model.integer
    @assert length(integers) == length(x0)
    nt1 = options.subsolver_options.nt
    subsolver_options = map(keys(nt1)) do k
        string(k) => nt1[k]
    end
    nl_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, Dict(subsolver_options)...)
    nt2 = options.nt
    solver_options = map(nt2) do k
        string(k) => nt2[k]
    end
    optimizer = JuMP.optimizer_with_attributes(
        Pavito.Optimizer,
        "cont_solver" => nl_solver, Dict(solver_options)...,
        "mip_solver" => JuMP.optimizer_with_attributes(Cbc.Optimizer),
    )
    problem, counter = get_jump_problem(
        model, x0; first_order = options.first_order,
        optimizer = optimizer, integers = integers,
    )
    return PavitoIpoptCbcWorkspace(model, problem, x0, integers, options, counter)
end
@params struct PavitoIpoptCbcResult
    minimizer
    minimum
    problem
    status
    fcalls::Int
end

function optimize!(workspace::PavitoIpoptCbcWorkspace)
    @unpack problem, options, counter = workspace
    counter[] = 0
    jump_problem = workspace.problem
    jump_model = jump_problem.model
    moi_model = jump_model.moi_backend
    MOI.optimize!(moi_model)
    minimizer = MOI.get(moi_model, MOI.VariablePrimal(), jump_problem.vars)
    objval = MOI.get(moi_model, MOI.ObjectiveValue())
    term_status = MOI.get(moi_model, MOI.TerminationStatus())
    primal_status = MOI.get(moi_model, MOI.PrimalStatus())
    return PavitoIpoptCbcResult(
        minimizer, objval, problem, (term_status, primal_status), counter[],
    )
end

struct PavitoIpoptCbcAlg{O} <: AbstractOptimizer
    options::O
end
PavitoIpoptCbcAlg(; kwargs...) = PavitoIpoptCbcAlg(kwargs)

function Workspace(model::VecModel, optimizer::PavitoIpoptCbcAlg, args...; kwargs...,)
    return PavitoIpoptCbcWorkspace(model, args...; kwargs...)
end
