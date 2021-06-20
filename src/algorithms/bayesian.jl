using AbstractGPs, IntervalArithmetic

mutable struct ZeroOrderGPSurrogate <: Function
    f::Function
    gp::GP
    X::Vector{Vector{Float64}}
    y::Union{Vector{Vector{Float64}}, Vector{Float64}}
    noise::Float64
    std_multiple::Float64
    mode::Symbol
    N::Int
end
function ZeroOrderGPSurrogate(
    f, x; kernel = SqExponentialKernel(),
    X = [x], y = [f(x)],
    noise = 1e-8, std_multiple = 3.0, mode = :interval,
)
    @assert noise > 0
    gp = GP(kernel)
    N = length(y[1])
    return ZeroOrderGPSurrogate(
        f, gp, X, y, noise, std_multiple, mode, N,
    )
end

function (s::ZeroOrderGPSurrogate)(x)
    if s.mode == :exact
        y = s.f(x)
        s.X = vcat(s.X, [x])
        s.y = vcat(s.y, [y])
  	    return Interval.(y, y)
    else
        if eltype(s.y) <: Real
            _m, _v = mean_and_var(posterior(
                s.gp(s.X, s.noise), s.y,
            ), [x])
            m, v = _m[1], _v[1]
        else
            _gp = s.gp(s.X, s.noise)
            ms_vs = map(1:s.N) do i
                mean_and_var(posterior(_gp, getindex.(s.y, i)), [x])
            end
            m = reduce(vcat, getindex.(ms_vs, 1))
            v = reduce(vcat, getindex.(ms_vs, 2))
        end
        r = s.std_multiple .* sqrt.(v)
        return Interval.(m .- r, m .+ r)
    end
end

function surrogate(f, x; kwargs...)
    s = ZeroOrderGPSurrogate(f, x; mode = :exact, kwargs...)
    s(x)
    s.mode = :interval
    return s
end

_lower_f(s) = x -> getproperty.(s(x), :lo)
_equality_f(s) = x -> begin
    t = s(x)
    return [getproperty.(t, :lo); .- getproperty(t, :hi)]
end

function surrogate_model(vecmodel::VecModel; kwargs...)
    surrogates = Function[]
    x0 = getinit(vecmodel)
    if :expensive in vecmodel.objective.flags
        s = surrogate(vecmodel.objective.f, copy(x0); kwargs...)
        push!(surrogates, s)
        obj = Objective(
            _lower_f(s),
            vecmodel.objective.multiple,
            vecmodel.objective.flags,
        )
    else
        obj = vecmodel.objective
    end
    expensive_eq_constraints = filter(vecmodel.eq_constraints.fs) do c
        :expensive in c.flags
    end
    cheap_eq_constraints = filter(vecmodel.eq_constraints.fs) do c
        @assert c isa EqConstraint
        !(:expensive in c.flags)
    end
    eq_constraints = VectorOfFunctions(cheap_eq_constraints)
    ineq_constraints1 = mapreduce(vcat, expensive_eq_constraints; init = Union{}[]) do c
        @assert c isa EqConstraint
        s = surrogate(c, copy(x0); kwargs...)
        push!(surrogates, s)
        return IneqConstraint(
            _equality_f(s), [zero.(c.rhs); zero.(c.rhs)], c.dim * 2, c.flags,
        )
    end
    ineq_constraints2 = map(vecmodel.ineq_constraints.fs) do c
        @assert c isa IneqConstraint
        if :expensive in c.flags
            s = surrogate(c, copy(x0); kwargs...)
            push!(surrogates, s)
            return IneqConstraint(
                _lower_f(s), zero.(c.rhs), c.dim, c.flags,
            )
        else
            return c
        end
    end
    ineq_constraints = VectorOfFunctions(
        vcat(ineq_constraints1, ineq_constraints2),
    )
    return VecModel(
        obj, eq_constraints, ineq_constraints, vecmodel.box_min,
        vecmodel.box_max, vecmodel.init, vecmodel.integer,
    ), surrogates
end

function update_surrogates!(model, surrogates, x)
    foreach(surrogates) do s
        s.mode = :exact
    end
    o = model.objective(x)
    i = model.ineq_constraints(x)
    e = model.eq_constraints(x)
    foreach(surrogates) do s
        s.mode = :interval
    end
    return o, i, e
end

@params struct BayesOptOptions
    sub_options
    maxiter::Int
    initialize::Bool
    ninit::Int
    ctol::Float64
    ftol::Float64
    postoptimize::Bool
    nt::NamedTuple
end
function BayesOptOptions(;
    sub_options = IpoptOptions(),
    initialize = true,
    ninit = 10,
    maxiter = 100,
    std_multiple = 3.0,
    kernel = SqExponentialKernel(),
    noise = 1e-8,
    ctol = 1e-4,
    ftol = 1e-4,
    postoptimize = false,
)
    return BayesOptOptions(
        sub_options,
        maxiter,
        initialize,
        ninit,
        ctol,
        ftol,
        postoptimize,
        (
            std_multiple = std_multiple,
            kernel = kernel,
            noise = noise,
        ),
    )
end

@params mutable struct BayesOptWorkspace <: Workspace
    model::VecModel
    sub_workspace
    x0::AbstractVector
    options::BayesOptOptions
    surrogates::AbstractVector
end

@params struct BayesOptResult <: AbstractResult
    minimizer
    minimum
    niters::Int
    sub_result::AbstractResult
    surrogates
end

# ScaledSobolIterator adapted from BayesianOptimization.jl
struct ScaledSobolIterator{T,D}
    lowerbounds::Vector{T}
    upperbounds::Vector{T}
    N::Int
    seq::SobolSeq{D}
end

"""
    ScaledSobolIterator(lowerbounds, upperbounds, N;
                        seq = SobolSeq(length(lowerbounds)))
Returns an iterator over `N` elements of a Sobol sequence between `lowerbounds`
and `upperbounds`. The first `N` elements of the Sobol sequence are skipped for
better uniformity (see https://github.com/stevengj/Sobol.jl)
"""
function ScaledSobolIterator(lowerbounds, upperbounds, N;
                             seq = SobolSeq(length(lowerbounds)))
    N > 0 && skip(seq, N)
    ScaledSobolIterator(lowerbounds, upperbounds, N, seq)
end
Base.length(it::ScaledSobolIterator) = it.N
@inline function Base.iterate(it::ScaledSobolIterator, s = 1)
    s == it.N + 1 && return nothing
    Sobol.next!(it.seq, it.lowerbounds, it.upperbounds), s + 1
end

function optimize!(workspace::BayesOptWorkspace)
    @unpack sub_workspace, surrogates = workspace
    @unpack options, x0 = workspace
    @unpack maxiter, ctol, ftol, postoptimize = options
    @unpack ninit, initialize = options
    smodel = sub_workspace.model
    minval = Inf
    minimizer = copy(x0)
    if initialize
        initializer = ScaledSobolIterator(getmin(smodel), getmax(smodel),
        ninit)
        x2 = x0
        for x1 in initializer
            sub_workspace.x0 .= x1
            r = optimize!(sub_workspace)
            if r.minimizer != x2
                ob, ineq, eq = update_surrogates!(smodel, surrogates, r.minimizer)
                feasible = all(ineq .<= ctol) && all(abs.(eq) .<= ctol)
                if feasible && ob <= minval
                    minval = ob
                    minimizer = copy(r.minimizer)
                end
            end
            x2 = r.minimizer
        end
    end
    sub_workspace.x0 .= minimizer
    r = optimize!(sub_workspace)
    m = r.minimum
    x = r.minimizer
    iter = 1
    while iter < maxiter
        iter += 1
        ob, ineq, eq = update_surrogates!(smodel, surrogates, x)
        feasible = all(ineq .<= ctol) && all(abs.(eq) .<= ctol)
        if feasible && ob <= minval
            minval = ob
            minimizer = copy(x)
        end
        converged = feasible && ob <= m + ftol
        converged && break
        r = optimize!(sub_workspace)
        m = r.minimum
        x = r.minimizer
    end
    ob, ineq, eq = update_surrogates!(smodel, surrogates, x)
    feasible = all(ineq .<= ctol) && all(abs.(eq) .<= ctol)
    if feasible && ob <= minval
        minval = ob
        minimizer = r.minimizer
    end
    if postoptimize
        foreach(surrogates) do s
            s.mode = :exact
        end
    end
    sub_workspace.x0 .= minimizer
    r = optimize!(sub_workspace)
    if postoptimize
        if r.minimum <= minval
            minval = r.minimum
            minimizer = copy(r.minimizer)
        end
        foreach(surrogates) do s
            s.mode = :interval
        end
    end
    return BayesOptResult(
        minimizer, minval, iter, r, surrogates,
    )
end

struct BayesOptAlg{A} <: AbstractOptimizer
    sub_alg::A
end

function Workspace(
    model::VecModel, optimizer::BayesOptAlg, x0::AbstractVector;
    options::BayesOptOptions, surrogates = nothing, kwargs...,
)
    if surrogates !== nothing
        smodel = model
    else
        smodel, surrogates = surrogate_model(model; options.nt...)
    end
    sub_workspace = Workspace(
        smodel,
        optimizer.sub_alg, x0;
        options = options.sub_options, kwargs...,
    )
    return BayesOptWorkspace(
        model,
        sub_workspace,
        x0,
        options,
        surrogates,
    )
end
