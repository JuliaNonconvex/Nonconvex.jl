"""
    MsOAlg

 A struct to wrap the algorithms from `MultistartOptimzation.jl`.
 - `alg`: the global optimization algorithm
"""
struct MsOAlg
    alg
end

"""
    MsOOptions

 A struct to store the options for `multistart_minimization` from `MultistartOptimzation.jl`.
 - `nt`: A named tuple to store all the options
"""
@params struct MsOOptions
    nt::NamedTuple
end
function MsOOptions(;
    use_threads = true,
)
    MsOOptions(
        (;
            use_threads = use_threads
        )
    )
end

"""
    MsOWorkspace

A struct that stores all the necessary memory allocations needed for global optimization through `multistart_minimization` from `MultistartOptimzation.jl`. The following are the fields of the struct:
 - `model`: the original model to be minimized
 - `problem`: an instance of [`MultistartOptimization.MinimizationProblem`] that defines the minimization problem
 - `optimizer`: an instance of [`MsOAlg`](@ref) that stores the global optimization algorithm to be used
 - `options`: an instance of [`MsOOptions`](@ref) that resembles options for `multistart_minimization`
 - `local_method`: the local optimization algorithm used
"""
@params mutable struct MsOWorkspace <: Workspace
    model::Model
    problem::MultistartOptimization.MinimizationProblem
    optimizer::MsOAlg
    options::MsOOptions
    local_method
end
function MsOWorkspace(
    model::Model, optimizer::MsOAlg,
    local_method, options = MsOOptions(),
)
    problem = getMsO_problem(optimizer.alg, model, local_method)
    
    return MsOWorkspace(model, problem, optimizer, options, local_method)
end


function getMsO_problem(alg, model::Model, local_method)
    obj = model.objective.f

    xlb = getmin(model)
    xub = getmax(model)

    if length(xlb) == 0 || length(xub) == 0
        error("Must specify lower and upper bounds for objective.")
    end

    problem = MultistartOptimization.MinimizationProblem(obj, xlb, xub)

    return problem
end

"""
    MsOResult

 A summary result struct returned by [`optimize_MsO`](@ref) that includes the following fields:
 - `minimizer`: the best solution found
 - `minimum`: the optimal value
 - `problem`: the minimization problem
 - `optimizer`: the global optimization algorithm
 - `local_method`: the local optimization algorithm used
"""
@params struct MsOResult
    minimizer
    minimum
    problem
    optimizer
    local_method
end


function optimize_MsO(workspace::MsOWorkspace)
    @unpack problem, optimizer, options, local_method = workspace

    result  = MultistartOptimization.multistart_minimization(
                        optimizer.alg, local_method, problem,
                        use_threads = options.nt.use_threads
                    )

    return MsOResult(
        result.location, result.value, problem, optimizer, local_method
    )
end


"""
get_starting_point(
    model::Model, optimizer::MsOAlg;
    local_method, local_maxiters = nothing,
    options = MsOOptions(), kwargs...,
)

Returns a starting point after performing a global search on the `objective` of the `model`.
"""
function get_starting_point(
    model::Model, optimizer::MsOAlg;
    local_method, local_maxiters = nothing,
    options = MsOOptions(), kwargs...,
    )

    if !isnothing(local_maxiters) && local_maxiters <= 0.0
        error("The number of maximum iterations in local search (local_maxiters) needs to be positive.")
    elseif !isnothing(local_maxiters)
        local_maxiters = convert(Int, local_maxiters)
    end

    if !isnothing(local_maxiters)
        local_method = MultistartOptimization.NLoptLocalMethod(local_method, maxeval = local_maxiters)
    else
        local_method = MultistartOptimization.NLoptLocalMethod(local_method)
    end

    mso_workspace = MsOWorkspace(model, optimizer, local_method, options)

    mso_result = optimize_MsO(mso_workspace)

    return mso_result
end
