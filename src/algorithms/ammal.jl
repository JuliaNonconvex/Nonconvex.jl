"""
    MMALagOptions

A struct that stores all the options of the MMA algorithms. Th following are the fields of `MMALagOptions`:
 - `maxiter`: the maximum number of inner iterations. For `MMA87`, there is 1 inner iteration per outer iteration.
 - `outer_maxiter`: the maximum number of outer iterations.
 - `maxinner`: the maximum number of inner iterations per outer iteration of [`MMA02`](@ref). Not applicable for [`MMA87`](@ref).
 - `tol`: a tolerance struct of type [`Tolerance`](@ref).
 - `s_init`: defined in the original [`MMA02`](@ref) paper.
 - `s_incr`: defined in the original [`MMA02`](@ref) paper.
 - `s_decr`: defined in the original [`MMA02`](@ref) paper.
 - `store_trace`: if true, a trace will be stored.
 - `auto_scale`: if set to `true`, the objective will be auto-scaled to make the norm of the objective and its gradient similar to the constriants'.
 - `dual_options`: the options passed to the dual optimizer from [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl).
 - `tonormalize`: if set to `true`, the dual optimizer will normalize the dual variables and perform gradient descent on the hyper-sphere.
 - `norm`: the norm type used to normalize the dual variables if `tonormalize` is set to `true`.
 - `callback_steps`: the [`MMALagCallback`](@ref) is called every `callback_steps` outer iterations. The callback increments the quadratic penalty factor and the linear weights are updated.
 - `quad_incr`: the factor by which the quadratic penalty is incremented every time the [`MMALagCallback`](@ref) is called.
 - `dual_step`: the dual solution will be incremented according to a gradient descent scheme which is `Descent(dual_step * nconstr)` by default where `nconstr` is the number of constraints in the problem.
 - `dual_opt`: the dual gradient descent scheme that the dual solution is incremented by. This defaults to `Descent(dual_step * nconstr)` where `dual_step` is the keyword argument passed in and `nconstr` is the number of constraints. Other examples include [`Nesterov`](@ref), [`Momentum`](@ref) and [`ExpDecay`](@ref).
 - `verbose`: when set to `true`, more information is printed to the screen.
"""
@with_kw mutable struct MMALagOptions{T, Ttol <: Tolerance, TSubOptions <: Optim.Options}
    maxiter::Int = 1000
    outer_maxiter::Int = 10^8
    maxinner::Int = 10
    tol::Ttol = Tolerance()
    s_init::T = 0.1
    s_incr::T = 1.2
    s_decr::T = 0.7
    store_trace::Bool = false
    show_trace::Bool = false
    auto_scale::Bool = false
    keep_best::Bool = false
    dual_options::TSubOptions = Optim.Options(allow_f_increases = false, iterations = 1000, outer_iterations=1000)
    tonormalize::Bool = true
    norm = 2
    callback_steps::Int = 3
    quad_incr::T = 1.2
    dual_step::T = 0.01
    dual_opt = nothing
    verbose::Bool = false
end

@params struct MMALag <: AbstractOptimizer
    parent::AbstractOptimizer
end
MMALag() = MMALag(MMA02())

default_options(model::MMALagModel, alg::MMA87) = MMALagOptions()
default_options(model::MMALagModel, alg::MMA02) = MMALagOptions()
default_options(model::VecModel, alg::MMALag) = MMALagOptions()

@params mutable struct MMALagCallback <: Function
    previter::Int
    prevf::Real
    stochopt
end
function MMALagCallback(options, nconstr::Int=1)
    @unpack dual_opt, dual_step = options
    #opt = ExpDecay(0.1, 0.9, 1, 1e-6)
    #opt = Descent(nconstr * 0.01)
    nconstr = 1
    if dual_opt === nothing
        dual_opt = Descent(nconstr * dual_step)
    end
    return MMALagCallback(0, Inf, dual_opt)
end

function (callback::MMALagCallback)(workspace::MMAWorkspace)
    @unpack dualmodel, optimizer, options, solution = workspace
    @unpack convstate = solution
    @unpack tol = options
    infeastol = tol.infeas
    model = dualmodel |> getparent |> getparent
    callbackrun = workspace.outer_iter == 2 ||
        workspace.iter >= options.maxiter ||
        workspace.outer_iter >= options.outer_maxiter ||
        convstate.converged && convstate.infeas > infeastol ||
        workspace.outer_iter - callback.previter > options.callback_steps

    if callbackrun
        mmalagupdate!(
            model,
            workspace,
            solution.prevf,
            callback.stochopt,
            first = workspace.outer_iter == 2,
        )
        callback.previter = workspace.outer_iter
    end
    return callbackrun
end

mmalagupdate!(model::VecModel, workspace, prevf, stochopt; first = false) = model
function mmalagupdate!(model::MMALagModel, workspace, prevf, stochopt; first = false)
    @unpack dualmodel, solution, optimizer, options, convcriteria = workspace
    @unpack convstate = solution
    @unpack tol, tonormalize = options
    normn = options.norm
    ftol, infeastol = tol.f, tol.infeas
    model = dualmodel |> getparent |> getparent

    f = solution.f
    g = getorigconstrval(model)

    N = length(getlinweights(model)) ÷ 2
    if tonormalize
        dims = getdim.(getineqconstraints(getparent(model)).fs)
        nc = length(dims)

        λ = getlinweights(model)[1:N]
        start = 1
        for i in 1:nc
            inds = start:start+dims[i]-1
            λ[inds] .= getlinweights(model)[inds] .* solution.λ[i] .+ getlinweights(model)[inds .+ N] .* solution.λ[i .+ nc]
            start += dims[i]
        end
        normalize!(λ, normn)

        direction = copy(g)
        start = 1
        for i in 1:nc
            inds = start:start+dims[i]-1
            normval = norm(direction[inds], normn)
            if normval != 0
                direction[inds] .= direction[inds] ./ normval
            end
            start += dims[i]
        end
        apply!(stochopt, λ, direction)

        newλ = similar(λ)
        start = 1
        for i in 1:nc
            inds = start:start+dims[i]-1
            if length(inds) == 1
                newλ[start] = 1
            else
                temp = max.(0, (λ[inds] .+ direction[inds]))
                newλ[inds] = temp
                normval = norm(temp, normn)
                if normval != 0
                    newλ[inds] .= newλ[inds] ./ normval
                end
            end
            start += dims[i]
        end
        setlinweights!(model, [λ; newλ])
    else
        dims = getdim.(getineqconstraints(getparent(model)).fs)
        nc = length(dims)

        λ = getlinweights(model)[1:N]
        start = 1
        for i in 1:nc
            inds = start:start+dims[i]-1
            λ[inds] .= getlinweights(model)[inds] .* solution.λ[i] .+ getlinweights(model)[inds .+ N] .* solution.λ[i + nc]
            start += dims[i]
        end
        direction = copy(g)#.* (1 .+ rand.())
        apply!(stochopt, λ, direction)
        newλ = max.(0, (λ .+ direction))
        setlinweights!(model, [λ; newλ])
    end

    c = getquadweight(model)
    if convstate.infeas <= infeastol
        newc = c
    else
        newc = min(options.quad_incr * c, 1e6)
    end
    setquadweight!(model, newc)
    if debugging[]
        @show getquadweight(model)

        println()
        #@show maximum(g)
        #@show convstate.infeas
        #@show getoptimalx(dualmodel), solution.f
        #@show newλ, newc
    end

    approxmodel = getparent(dualmodel)
    optimalx = getoptimalx(dualmodel)
    aggfg, agg∇fg = value_jacobian(
        approxmodel.objective_ineq_constraints,
        optimalx,
    )
    workspace.fcalls += 1

    # Update the convex approximation at the current x
    updateapprox!(dualmodel, optimalx, aggfg, agg∇fg)
    
    # Update the objective and constraint values and gradients in `solution`
    updatefg!(solution, aggfg, agg∇fg)

    correctsolution!(solution, model, options)

    return model
end

function init!(model::MMALagModel)
    set_objective_multiple!(model, 1)
    setquadweight!(model, 1e-5)
    setlinweights!(model, 1)
    return model
end

function scalequadweight!(model::MMALagModel, s::Real)
    setquadweight!(model, getquadweight(model) * s)
    return model
end

function correctsolution!(solution::Solution, model::MMALagModel, options::MMALagOptions)
    @unpack convstate = solution
    @unpack tol = options
    infeastol = tol.infeas
    f = getorigobjval(model)
    if debugging[]
        @show f
    end
    g = getorigconstrval(model)
    infeas = max(maximum(g), 0)
    if debugging[]
        @show infeas
    end
    convstate.infeas = infeas
    convstate.converged = convstate.converged && infeas <= infeastol
    solution.f = f
    return solution
end

function Workspace(model::VecModel, alg::MMALag, x0::AbstractVector, args...; kwargs...)
    return MMAWorkspace(MMALagModel(model; kwargs...), alg.parent, x0, args...; kwargs...)
end
