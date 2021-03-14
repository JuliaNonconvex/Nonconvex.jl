@params struct AugLag <: AbstractOptimizer
    primaloptimizer::AbstractOptimizer
    dualoptimizer::AbstractOptimizer
end
function AugLag(;
    primaloptimizer = Optim.ConjugateGradient(linesearch=Optim.LineSearches.BackTracking(iterations = 10)),
    dualoptimizer = Optim.GradientDescent(linesearch=Optim.LineSearches.BackTracking(iterations = 10)),
)
    return AugLag(primaloptimizer, dualoptimizer)
end

@params struct AugLagOptions
    primaloptions
    dualoptions
    maxiter
    tol
    quadfactor
end
function AugLagOptions(alg::AugLag;
    primaloptions = alg.primaloptimizer isa MMA02 || alg.primaloptimizer isa MMA87 ? MMAOptions(maxiter = 100, tol = Tolerance(kkt = 1e-4)) : Optim.Options(outer_iterations = 10, iterations = 10),
    dualoptions = Optim.Options(outer_iterations = 10, iterations = 10),
    maxiter = 10,
    tol = Tolerance(),
    quadfactor = 10,
)
    return AugLagOptions(primaloptions, dualoptions, maxiter, tol, quadfactor)
end

function Solution(lagmodel::AugLagModel)
    prevx = copy(getmin(lagmodel))
    x = copy(prevx)
    prevf = Inf
    f = Inf
    λ = copy(getlinweights(lagmodel))
    g = copy(λ)
    convstate = ConvergenceState()
    return Solution(prevx, x, getlinweights(lagmodel), prevf, f, nothing, g, nothing, convstate)
end

@params mutable struct AugLagWorkspace <: Workspace
    model::Model
    lagmodel::AugLagModel
    x0::AbstractVector
    optimizer::AugLag
    options::AugLagOptions

    solution::Solution
    convcriteria::ConvergenceCriteria
    callback::Function

    trace::Trace
    outer_iter::Int
    iter::Int
	fcalls::Int
end
function AugLagWorkspace(
    model::AbstractModel,
    optimizer::AugLag,
    x0::AbstractVector;
    options::AugLagOptions = AugLagOptions(optimizer),
    convcriteria::ConvergenceCriteria = KKTCriteria(),
    plot_trace::Bool = false,
    show_plot::Bool = plot_trace,
    save_plot = nothing,
    callback::Function = plot_trace ? LazyPlottingCallback(; show_plot = show_plot, save_plot = save_plot) : NoCallback(),
    kwargs...,
)
    T = eltype(x0)
    lagmodel = AugLagModel(model; kwargs...)

    # Convergence
    solution = Solution(lagmodel)
    #assess_convergence!(solution, model, options.tol, convcriteria)

    # Trace
    trace = Trace([])

    # Iteraton counter
    fcalls, outer_iter, iter = 1, 0, 0

    return AugLagWorkspace(
        model,
        lagmodel,
        x0,
        optimizer,
        options,

        solution,
        convcriteria,
        callback,

        trace,
        outer_iter,
        iter,
        fcalls,
    )
end

Workspace(model::AbstractModel, alg::AugLag, x0::AbstractVector; kwargs...) = AugLagWorkspace(model, alg, x0; kwargs...)

function optimize!(workspace::AugLagWorkspace)
    @unpack lagmodel, solution, options, convcriteria = workspace
    @unpack callback, optimizer, trace = workspace
    @unpack x0, outer_iter, iter, fcalls = workspace
    @unpack primaloptimizer, dualoptimizer = optimizer
    @unpack maxiter, primaloptions, dualoptions, quadfactor = options
    @unpack prevx, x, λ = solution

    xl, xu = getmin(lagmodel), getmax(lagmodel)

    T = eltype(x)
    ni, nj = length(λ), length(x)

    # Original model
    model = lagmodel |> getparent

    auglag = getobjective(lagmodel)

    cb = (tr; kwargs...) -> begin
        solution = deepcopy(solution)
        solution.prevx .= solution.x
        solution.x .= getx(lagmodel)
        solution.λ .= getλ(lagmodel)
        solution.prevf = solution.f
        solution.f = getorigobjval(lagmodel)
        solution.g .= getorigconstrval(lagmodel)
        assess_convergence!(solution, lagmodel, options.tol, convcriteria)
        callback(solution)
        return hasconverged(solution)
	end

    primaloptimizerfunction(λ) = begin
        setlinweights!(auglag, λ)
        # Primal objective to be minimized
        # Calculates the objective and its gradient
        # Optim-compatible objective
        if primaloptimizer isa MMA02 || primaloptimizer isa MMA87
            primalmodel = Model()
            addvar!(primalmodel, xl, xu)
            set_objective!(primalmodel, getprimalobjective(auglag))
            result = optimize(primalmodel, primaloptimizer, clamp.(getx(lagmodel), xl .+ 1e-3, xu .- 1e-3), options = primaloptions, callback = cb, convcriteria = KKTCriteria())
            fcalls += result.fcalls
        else
            primalobj = getoptimobj(getprimalobjective(auglag), true)
            result = Optim.optimize(
                Optim.only_fg!(primalobj),
                xl,
                xu,
                clamp.(getx(lagmodel), xl .+ 1e-3, xu .- 1e-3),
                Optim.Fminbox(primaloptimizer),
                @set(primaloptions.callback = cb),
            )
            fcalls += result.f_calls
        end
        if getx(lagmodel) != result.minimizer
            primalobj(result.minimizer)
            setx!(lagmodel, result.minimizer)
        end
        if debugging[]
            @show getx(lagmodel)
        end
        return getx(lagmodel), getorigobjval(lagmodel), getorigconstrval(lagmodel)
    end

    # Dual objective to be minimized - original objective will be maximized
    # Calculates negative the objective and its gradient
    # Optim-compatible objective
    dualobj = getoptimobj(getdualobjective(auglag, primaloptimizerfunction), false)

    # Lower and upper bounds on the dual variables
    λl = zeros(ni) .+ 1e-10
    λu = fill(Inf, ni)

    # Solve the dual problem by minimizing negative the dual objective value
    for i in 1:maxiter
        setquadweight!(lagmodel, min(getquadweight(lagmodel) * quadfactor, 1e10))
        if debugging[]
            @show getquadweight(lagmodel)
        end
        λresult = Optim.optimize(
            Optim.only_fg!(dualobj),
            λl,
            λu,
            max.(copy(λ), 1e-3),
            Optim.Fminbox(dualoptimizer),
            dualoptions, #@set(dualoptions.callback = cb),
        )
        λ .= λresult.minimizer
        if debugging[]
            @show λ
            @show getx(lagmodel)
        end
        if λ != getλ(lagmodel)
            dualobj(λ)
            setλ!(lagmodel, λ)
        end
    end
    if debugging[]
        #@show getx(lagmodel)
    end

    @pack! workspace = iter, fcalls
    solution.x .= getx(lagmodel)
    solution.λ .= getλ(lagmodel)
    solution.f = getorigobjval(lagmodel)
    solution.g .= getorigconstrval(lagmodel)
    callback(solution, update = true)

    results = GenericResult(
        optimizer,
        x0,
        solution.x,
        solution.f,
        iter,
        iter == options.maxiter,
        options.tol,
        solution.convstate,
        fcalls,
    )
    return results
end
