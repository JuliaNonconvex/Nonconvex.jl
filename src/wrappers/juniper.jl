import .Juniper

@params struct JuniperIpoptOptions
    nt::NamedTuple
    subsolver_options
    first_order::Bool
end
function JuniperIpoptOptions(;
    first_order = true,
    linear_constraints = false,
    subsolver_options = IpoptOptions(
        first_order = first_order,
        linear_constraints = linear_constraints,
    ),
    kwargs...,
)
    first_order = !hasproperty(subsolver_options.nt, :hessian_approximation) ||
        subsolver_options.nt.hessian_approximation == "limited-memory"
    return JuniperIpoptOptions((; kwargs...), subsolver_options, first_order)
end

@params mutable struct JuniperIpoptWorkspace <: Workspace
    model::VecModel
    problem::JuMPProblem
    x0::AbstractVector
    integers::AbstractVector{<:Bool}
    options::JuniperIpoptOptions
    counter::Base.RefValue{Int}
end
function JuniperIpoptWorkspace(
    model::VecModel, x0::AbstractVector = getinit(model);
    options = JuniperIpoptOptions(), kwargs...,
)
    integers = model.integer
    @assert length(integers) == length(x0)
    nt1 = options.subsolver_options.nt
    subsolver_options = map(keys(nt1)) do k
        string(k) => nt1[k]
    end
    nl_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, Dict(subsolver_options)...)
    nt2 = options.nt
    solver_options = map(keys(nt2)) do k
        string(k) => nt2[k]
    end
    optimizer = JuMP.optimizer_with_attributes(
        Juniper.Optimizer, "nl_solver" => nl_solver, Dict(solver_options)...,
    )
    problem, counter = get_jump_problem(
        model, x0; first_order = options.first_order,
        optimizer = optimizer, integers = integers,
    )
    return JuniperIpoptWorkspace(model, problem, x0, integers, options, counter)
end
@params struct JuniperIpoptResult <: AbstractResult
    minimizer
    minimum
    problem
    status
    fcalls::Int
end

function optimize!(workspace::JuniperIpoptWorkspace)
    @unpack problem, options, x0, counter = workspace
    counter[] = 0
    jump_problem = workspace.problem
    jump_model = jump_problem.model
    moi_model = jump_model.moi_backend
    MOI.set.(Ref(moi_model), Ref(MOI.VariablePrimalStart()), workspace.problem.vars, x0)
    MOI.optimize!(moi_model)
    minimizer = MOI.get(moi_model, MOI.VariablePrimal(), jump_problem.vars)
    objval = MOI.get(moi_model, MOI.ObjectiveValue())
    term_status = MOI.get(moi_model, MOI.TerminationStatus())
    primal_status = MOI.get(moi_model, MOI.PrimalStatus())
    return JuniperIpoptResult(
        minimizer, objval, problem, (term_status, primal_status), counter[],
    )
end

struct JuniperIpoptAlg{O} <: AbstractOptimizer
    options::O
end
JuniperIpoptAlg(; kwargs...) = JuniperIpoptAlg(kwargs)

function Workspace(model::VecModel, optimizer::JuniperIpoptAlg, args...; kwargs...,)
    return JuniperIpoptWorkspace(model, args...; kwargs...)
end
