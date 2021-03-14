# MMA Algorithms

"""
    MMA87

The original method of moving asymptotes (MMA) algorithm from the [1987 paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207).
"""
@params struct MMA87 <: AbstractOptimizer
    dualoptimizer::AbstractOptimizer
end
function MMA87(; dualoptimizer = Optim.GradientDescent())
    return MMA87(dualoptimizer)
end

"""
    MMA02

The globally convergent method of moving asymptotes (MMA) algorithm from the [2002 paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822).
"""
@params struct MMA02 <: AbstractOptimizer
    dualoptimizer::AbstractOptimizer
end
function MMA02(; dualoptimizer = Optim.GradientDescent())
    return MMA02(dualoptimizer)
end

"""
    Trace

A struct that stores the history of solutions.
"""
struct Trace
    tr::Vector
end

"""
    GenericResult

A summary result struct returned by [`optimize`](@ref) that includes the following fields:
 - `optimizer`: the optimization algorithm used
 - `initial_x`: the initial primal solution
 - `minimizer`: the best solution found
 - `minimum`: the optimal value
 - `iter`: number of inner iterations run
 - `maxiter_reached`: true if the algorithm stopped due to reaching the maximum number of iterations
 - `tol`: an instance of [`Tolerance`](@ref) that specifies the convergence tolerances
 - `convstate`: an instance of [`ConvergenceCriteria`](@ref) that summarizes the convergenc state of the best solution found
 - `fcalls`: the number of times the objective and constraint functions were called during the optimization
"""
@params mutable struct GenericResult{T}
    optimizer
    initial_x::AbstractVector{T}
    minimizer::AbstractVector{T}
    minimum::T
    iter::Int
    maxiter_reached::Bool
    tol::Tolerance
    convstate::Union{ConvergenceState, Nothing}
    fcalls::Int
end

"""
    MMAWorkspace

A struct that stores all the intermediate states and memory allocations needed for the optimization. The following are the fields of the struct:
 - `model`: the original instance of [`MMAModel`](@ref) being optimized.
 - `dualmodel`: an instance of [`DualMMAModel`](@ref) which resembles the dual of the MMA approximation of the original model.
 - `x0`: the initial primal solution.
 - `solution`: an instance of [`Solution`](@ref) that resembles the current solution in the optimization algorithm.
 - `σ`: the workspace move limit, explained in [`MMAApprox`](@ref).
 - `ρ`: the `ρ` parameter as explained in [`MMAApprox`](@ref).
 - `tempx`: a temporary vector used to store the 2nd previous primal solution.
 - `options`: an instance of [`MMAOptions`](@ref) that resembles the options of the MMA algorithm.
 - `convcriteria`: an instance of [`ConvergenceCriteria`](@ref) that specifies the convergence criteria of the MMA algorithm.
 - `callback`: a function that is called on `solution` in every iteration of the algorithm. This can be used to store information about the optimization process.
 - `optimizer`: an instance of [`AbstractOptimizer`](@ref) such as `MMA87()` or `MMA02()` that specifies the variant of MMA used to optimize the model.
 - `suboptimizer`: the dual optimization algorithm used to optimize the barrier problem. This should be an [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl) optimizer.
 - `tracing`: a Boolean that when true stores the trace of solutions evaluated during the optimization process.
 - `outer_iter`: the current outer iteration index. [`MMA87`](@ref) has 1 inner iteration per outer iteration. [`MMA02`](@ref) does multiple inner iterations per outer iteration.
 - `iter`: the current inner iteration index.
 - `fcalls`: the current number of objective and constraint function calls.
"""
@params mutable struct MMAWorkspace <: Workspace
    model::AbstractModel
    dualmodel::MMADualModel
    x0::AbstractVector
    solution::Solution
    σ::AbstractVector
    ρ::AbstractVector
    tempx::AbstractVector
    convcriteria::ConvergenceCriteria
    callback::Function
    optimizer::AbstractOptimizer
    options
    trace::Trace
    outer_iter::Int
    iter::Int
	fcalls::Int
end
function MMAWorkspace(
    model::AbstractModel,
    optimizer::AbstractOptimizer,
    x0::AbstractVector;
    options = default_options(model, optimizer),
    convcriteria::ConvergenceCriteria = KKTCriteria(),
    plot_trace::Bool = false,
    show_plot::Bool = plot_trace,
    save_plot = nothing,
    callback::Function = plot_trace ? LazyPlottingCallback(; show_plot = show_plot, save_plot = save_plot) : NoCallback(),
    kwargs...,
)
    T = eltype(x0)
    init!(model)
    dualmodel = MMADualModel(MMAApproxModel(model, x0; kwargs...))

    # Convergence
    λ = ones(getdim(getineqconstraints(model)))
    solution = Solution(dualmodel, λ)
    assess_convergence!(solution, model, options.tol, convcriteria)
    correctsolution!(solution, model, options)

    # Trace
    trace = Trace([])

    # Iteraton counter
    fcalls, outer_iter, iter = 1, 0, 0
    σ = similar(x0)
    ρ = similar(λ, length(λ) + 1)
    tempx = similar(x0)

    return MMAWorkspace(
        model,
        dualmodel,
        x0,
        solution,
        σ,
        ρ,
        tempx,
        convcriteria,
        callback,
        optimizer,
        options,
        trace,
        outer_iter,
        iter,
        fcalls,
    )
end

default_options(model::Model, alg::MMA87) = MMAOptions()
default_options(model::Model, alg::MMA02) = MMAOptions()

init!(model::Model) = model

function Workspace(model::AbstractModel, optimizer::Union{MMA87, MMA02}, args...; kwargs...)
    return MMAWorkspace(model, optimizer, args...; kwargs...)
end

"""
```
optimize(
    model::AbstractModel,
    optimizer::AbstractOptimizer = MMA02(),
    x0::AbstractVector;
    options,
    convcriteria::ConvergenceCriteria = KKTCriteria(),
    plot_trace::Bool = false,
    callback::Function = plot_trace ? LazyPlottingCallback() : NoCallback(),
    kwargs...,
)
```

Optimizes `model` using the algorithm `optimizer`, e.g. an instance of [`MMA87`](@ref) or [`MMA02`](@ref). `x0` is the initial solution. The keyword arguments are:
 - `options`: used to set the optimization options. It is an instance of [`MMAOptions`](@ref) for [`MMA87`](@ref) and [`MMA02`](@ref).
 - `convcriteria`: an instance of [`ConvergenceCriteria`](@ref) that specifies the convergence criteria of the MMA algorithm.
 - `plot_trace`: a Boolean that if true specifies the callback to be an instance of [`PlottingCallback`](@ref) and plots a live trace of the last 50 solutions.
 - `callback`: a function that is called on `solution` in every iteration of the algorithm. This can be used to store information about the optimization process.

 The details of the MMA optimization algorithms can be found in the original [1987 MMA paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207) and the [2002 paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822).
"""
function optimize(args...; kwargs...)
    workspace = Workspace(args...; kwargs...)
    return optimize!(workspace)
end

function optimize!(workspace::MMAWorkspace)
    @unpack dualmodel, solution, convcriteria = workspace
    @unpack callback, optimizer, options, trace = workspace
    @unpack x0, σ, ρ, outer_iter, iter, fcalls = workspace
    @unpack dualoptimizer = optimizer
    @unpack dual_options, maxiter, outer_maxiter, auto_scale = options 
    @unpack prevx, x, g, λ = solution
    best_solution = deepcopy(solution)

    T = eltype(x)
    ni, nj = length(λ), length(x)

    # Original model
    model = dualmodel |> getparent |> getparent

    # Approximate MMA model
    approxmodel = getparent(dualmodel)

    # Dual objective to be maximized
    dualobj = getobjective(dualmodel)

    # Optim-compatible objective
    # Calculates negative the dual objective and its gradient
    optimobj = getoptimobj(dualobj, false)

    # Lower and upper bounds on the dual variables
    λl = zeros(ni) .+ 1e-10
    λu = fill(Inf, ni)
    λ .= 1

    # Callback, e.g. a trace plotting callback
    callback(solution)

    # Initialize the workspace's trust region σ
    initializeσ!(workspace)
    # Initialize the workspace's lift ρ, only used in MMA02
    initializeρ!(workspace)

    while (!hasconverged(solution) || outer_iter == 0) && iter < maxiter && outer_iter < outer_maxiter
        outer_iter += 1

        # Adapt workspace's trust region σ using s_incr and s_decr
        outer_iter > 2 && updateσ!(workspace)
        # Decrease the workspace's lift ρ, only used in MMA02
        decreaseρ!(workspace, optimizer)

        # Update the model's σ
        setσ!(dualmodel, σ)
        # Update the model's ρ, only used in MMA02
        setρ!(dualmodel, ρ)

        # Update the convex approximation at the current x
        updateapprox!(dualmodel)

        # Assume the convex approximation is not an upper bound approximation at x
        upperbounded = false
        prev_iter = iter
        while (!hasconverged(solution) || outer_iter == 1) && !upperbounded && iter < maxiter && iter < prev_iter + options.maxinner
            iter += 1
            # Solve the dual problem by minimizing negative the dual objective value
            if length(λ) > 0
                λ .= 1.0
                λ .= Optim.optimize(
                    Optim.only_fg!(optimobj),
                    λl,
                    λu,
                    λ,
                    Optim.Fminbox(dualoptimizer),
                    dual_options,
                ).minimizer
            end
            if debugging[]
                @show λ
            end

            # Revaluate the optimal primal solution at the optimal dual solution
            dualobj(λ)

            # Update the x vector to the new one and update the previous 2 xs
            updatex!(workspace)

            # Evaluate the approximate objective and constraints at the optimal x
            approxfg = getapproxfg(dualobj)

            # Evaluates the exact objective and constraints and their gradients at the optimal x
            optimalx = getoptimalx(dualmodel)
            fg, ∇fg = value_jacobian(
                approxmodel.objective_ineq_constraints,
                optimalx,
            )
            fcalls += 1

            # Scale the objective appropriately
            # Reduces the chance of overflow or underflow in λs 
            # Helps with the convergence
            auto_scale && iter < 5 && scaleobjective!(model, fg, ∇fg, approxfg, λ)
            #iter < 10 && scaleobjective!(model, fg, ∇fg, approxfg, λ)

            # Check if the approximation is an upper bound at the current x.
            # If not, increase ρ. Only used in MMA02
            if iter > 1
                upperbounded = increaseρ!(workspace, fg, approxfg, optimizer)
                upperbounded || setρ!(dualmodel, ρ)
            end
            # Update the convex approximation at the current x
            updateapprox!(dualmodel, optimalx, fg, ∇fg)
            # Update the objective and constraint values and gradients in `solution`
            updatefg!(solution, fg, ∇fg)

            # Check if the algorithm has converged
            assess_convergence!(solution, model, options.tol, convcriteria)

            # Callback, e.g. a trace plotting callback
            callback(solution)
        end
        @pack! workspace = outer_iter, iter, fcalls

        if options.keep_best
            best_solution = get_best_solution(solution, best_solution, options)
        else
            best_solution = deepcopy(solution) 
        end
        hasconverged(best_solution) && break

        # Print some trace if flag is on
        # @mmatrace()
    end
    @pack! workspace = outer_iter, iter, fcalls
    workspace.solution = best_solution

    # Reset the objective scaling factor to 1
    set_objective_multiple!(model, 1)
    callback(best_solution, update = true)
    
    results = GenericResult(
        optimizer,
        x0,
        best_solution.x,
        best_solution.f,
        iter,
        iter == options.maxiter,
        options.tol,
        best_solution.convstate,
        fcalls,
    )
    return results
end

function assess_convergence!(w::Workspace)
    assess_convergence!(w.solution, w.model, w.options.tol, w.convcriteria)
    correctsolution!(w.solution, w.model, w.options)
    return w
end

"""
    scaleobjective!(model::AbstractModel, fg, ∇fg, approxfg, λ)

Scales the objective of `model` using the ratio of the ∞ norms of the constraint and objective values and gradients. After scaling, the ∞ norms will be the same.
"""
function scaleobjective!(model::AbstractModel, fg, ∇fg, approxfg, λ)
    @views normratio = norm(∇fg[2:end,:], Inf) / norm(∇fg[1,:], Inf)
    @views normratio = max(normratio, norm(fg[2:end], Inf) / abs(fg[1]))
    normratio = isfinite(normratio) ? normratio : one(normratio)
    set_objective_multiple!(model, get_objective_multiple(model) * normratio)
    fg[1] = fg[1] * normratio
    approxfg[1] = approxfg[1] * normratio
    ∇fg[1,:] .= ∇fg[1,:] .* normratio
    scalequadweight!(model, normratio)
    λ .*= normratio
    return model
end

scalequadweight!(model::Model, s::Real) = model

function correctsolution!(solution::Solution, model::Model, options::MMAOptions)
    solution.f = solution.f / get_objective_multiple(model)
    return solution
end

function get_best_solution(solution::Solution, best_solution::Solution, options)
    best_infeas = best_solution.convstate.infeas
    if max(solution.convstate.infeas, solution.convstate.kkt_residual) <= 
        max(best_solution.convstate.infeas, best_solution.convstate.kkt_residual)
        best_solution = deepcopy(solution)
    end
    return best_solution
end

function updatex!(workspace::MMAWorkspace)
    @unpack tempx, solution, dualmodel = workspace
    @unpack prevx, x = solution
    tempx .= prevx
    prevx .= x
    x .= getoptimalx(dualmodel)
    return workspace
end

function updatefg!(solution::Solution, fg::AbstractVector, ∇fg::AbstractMatrix)
    solution.prevf = solution.f
    solution.f = fg[1]
    solution.g .= fg[2:end]
    solution.∇f .= ∇fg[1, :]
    solution.∇g .= ∇fg[2:end, :]
    return solution
end

function initializeσ!(workspace::MMAWorkspace)
    @unpack σ, optimizer, options, model = workspace
    T = eltype(σ)
    map!(σ, 1:length(σ)) do j
        diff = getmax(model, j) - getmin(model, j)
        if !isfinite(diff)
            diff = 1000 * one(T)
        end
        return options.s_init * diff
    end
    return workspace
end

function initializeρ!(workspace::MMAWorkspace)
    @unpack ρ = workspace
    ρ .= 0
    return workspace
end

function updateσ!(workspace::MMAWorkspace)
    @unpack σ, solution, tempx, optimizer, options, model = workspace
    @unpack x, prevx = solution
    T = eltype(x)
    map!(σ, 1:length(x)) do j
        if x[j] == prevx[j] || prevx[j] == tempx[j]
            d = σ[j]
        elseif xor(x[j] > prevx[j], prevx[j] > tempx[j])
            d = σ[j] * options.s_decr
        else
            d = σ[j] * options.s_incr
        end
        diff = getmax(model, j) - getmin(model, j)
        if !isfinite(diff)
            diff = 1000 * one(T)
        end
        min = diff/100
        max = 10diff
        if d <= min
            return min
        elseif d >= max
            return max
        else
            return d
        end
    end
    return workspace
end

decreaseρ!(workspace::MMAWorkspace, ::MMA87) = workspace
function decreaseρ!(workspace::MMAWorkspace, ::MMA02)
    @unpack ρ = workspace
    ρ .= max.(ρ ./ 10, 1e-5)
    return workspace
end

increaseρ!(workspace::MMAWorkspace, ::AbstractVector, ::AbstractVector, ::MMA87) = true
function increaseρ!(
    workspace::MMAWorkspace,
    fg::AbstractVector,
    approxfg::AbstractVector,
    ::MMA02,
)
    @unpack dualmodel, solution, σ, ρ = workspace
    @unpack x, prevx = solution
    T = eltype(x)
    w = zero(T)
    for j in 1:length(x)
        diff2 = (x[j] - prevx[j])^2
        w += diff2 / (σ[j]^2 - diff2)
    end
    w /= 2

    upperbounded = true
    for i in 1:length(ρ)
        if fg[i] > approxfg[i]
            upperbounded = false
            ρ[i] = min(10ρ[i], 1.1 * (ρ[i] + (fg[i] - approxfg[i]) / w))
        end
    end
    return upperbounded
end

function getoptimobj(obj, minimize = true)
    optimobj(z) = optimobj(1.0, nothing, z)
    function optimobj(F, G, z)
        if G !== nothing
            val, grad = value_gradient(obj, z)
            if minimize
                G[:] .= grad
            else
                G[:] .= .-grad
            end
            if F !== nothing
                if minimize
                    return val
                else
                    return -val
                end
            end
        end
        # No gradient necessary, just return the log joint.
        if F !== nothing
            if minimize
                return obj(z)
            else
                return -obj(z)
            end
        end
        return nothing
    end
    return optimobj
end
