using AbstractGPs, IntervalArithmetic

mutable struct ZeroOrderGPSurrogate <: Function
    f::Function
    gps::Vector{<:GP}
    X::Vector{Vector{Float64}}
    y::Union{Vector{Vector{Float64}}, Vector{Float64}}
    noise::Float64
    std_multiple::Float64
    mode::Symbol
    N::Int
    every::Int
    skip::Int
    last::Union{Symbol, Int}
    fit_prior::Bool
    fit_noise::Bool
    counter::Int
end
function ZeroOrderGPSurrogate(
    f, x; kernel = SqExponentialKernel(), X = [x], y = [f(x)], noise = 1e-8,
    std_multiple = 3.0, mode = :interval, every = 5, skip = 10, last = :all,
    fit_prior = true, fit_noise = false,
)
    @assert noise > 0
    gps = [GP(AbstractGPs.ConstMean(0.0), kernel) for _ in 1:length(y[1])]
    N = length(y[1])
    return ZeroOrderGPSurrogate(
        f, gps, X, y, noise, std_multiple, mode, N, every, skip, last,
        fit_prior, fit_noise, 0,
    )
end

function get_x_lb_ub(k::Kernel)
    x, un = flatten(k)
    lb = fill(-Inf, length(x))
    ub = fill(Inf, length(x))
    return x, lb, ub, un
end

function fit_mle!(s, gps, X, y, noise, last, fit_noise)
    xs_uns_1 = flatten.(getproperty.(gps, :mean))
    xs_lbs_ubs_uns_2 = get_x_lb_ub.(getproperty.(gps, :kernel))

    x1 = mapreduce(vcat, xs_uns_1) do x_un
        x_un[1]
    end
    x2 = mapreduce(vcat, xs_lbs_ubs_uns_2) do x_lb_ub_un
        x_lb_ub_un[1]
    end
    lb2 = mapreduce(vcat, xs_lbs_ubs_uns_2) do x_lb_ub_un
        x_lb_ub_un[2]
    end
    ub2 = mapreduce(vcat, xs_lbs_ubs_uns_2) do x_lb_ub_un
        x_lb_ub_un[3]
    end
    if fit_noise
        x = [x1; x2; noise]
        lb = [fill(-Inf, length(x1)); lb2; 1e-8]
        ub = [fill(Inf, length(x1)); ub2; Inf]
    else
        x = [x1; x2]
        lb = [fill(-Inf, length(x1)); lb2]
        ub = [fill(Inf, length(x1)); ub2]
    end
    if last == :all
        _X = X
        _y = y
    else
        _X = X[max(end-last+1, 1):end]
        _y = y[max(end-last+1, 1):end]
    end
    obj(θ) = begin
        if _y isa Vector{<:Vector}
            offset = 0
            return -sum(map(1:length(_y[1])) do i
                l1 = length(xs_uns_1[i][1])
                un1 = xs_uns_1[i][2]
                l2 = length(xs_lbs_ubs_uns_2[i][1])
                un2 = xs_lbs_ubs_uns_2[i][4]
                _gp = GP(un1(θ[offset+1:offset+l1]), un2(θ[offset+l1+1:offset+l1+l2]))
                offset += l1 + l2
                if fit_noise
                    return logpdf(_gp(_X, θ[end]), getindex.(_y, i))
                else
                    return logpdf(_gp(_X, noise), getindex.(_y, i))
                end
            end)
        else
            l1 = length(xs_uns_1[1][1])
            un1 = xs_uns_1[1][2]
            l2 = length(xs_lbs_ubs_uns_2[1][1])
            un2 = xs_lbs_ubs_uns_2[1][4]
            _gp = GP(un1(θ[1:l1]), un2(θ[l1+1:l1+l2]))
            if fit_noise
                return -logpdf(_gp(_X, θ[end]), _y)
            else
                return -logpdf(_gp(_X, noise), _y)
            end
        end
    end
    m = Model(obj)
    addvar!(m, lb, ub)
    options = IpoptOptions(max_iter = 20, print_level = 0)
    res = optimize(m, IpoptAlg(), x, options = options)
    x = res.minimizer
    offset = 0
    gps, noise = map(1:length(_y[1])) do i
        l1 = length(xs_uns_1[i][1])
        un1 = xs_uns_1[i][2]
        l2 = length(xs_lbs_ubs_uns_2[i][1])
        un2 = xs_lbs_ubs_uns_2[i][4]
        out = GP(un1(x[offset+1:offset+l1]), un2(x[offset+l1+1:offset+l1+l2]))
        offset += l1 + l2
        return out
    end, fit_noise ? x[end] : noise
    s.gps = gps
    s.noise = noise
    return s
end
function ChainRulesCore.rrule(::typeof(fit_mle!), s, gp, X, y, noise, last, fit_noise)
    return fit_mle!(s, gp, X, y, noise, last, fit_noise), _ -> begin
        return ntuple(_ -> NoTangent(), Val(8))
    end
end

function (s::ZeroOrderGPSurrogate)(x)
    if s.mode == :exact
        y = s.f(x)
        s.X = vcat(s.X, [x])
        s.y = vcat(s.y, [y])
        s.counter += 1
        if s.fit_prior && s.counter > s.skip * s.every && (s.counter % s.every) == 0
            fit_mle!(s, s.gps, s.X, s.y, s.noise, s.last, s.fit_noise)
        end
  	    return Interval.(y, y)
    else
        if eltype(s.y) <: Real
            _m, _v = mean_and_var(posterior(
                s.gps[1](s.X, s.noise), s.y,
            ), [x])
            m, v = _m[1], _v[1]
        else
            ms_vs = map(1:s.N) do i
                _gp = s.gps[i](s.X, s.noise)
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
    fit_prior = true,
    fit_noise = false,
    every = 2,
    skip = 2,
    last = :all,
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
            every = every,
            skip = skip,
            last = last,
            fit_prior = fit_prior,
            fit_noise = fit_noise,
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
    lb, ub = getmin(smodel), getmax(smodel)
    if initialize && !(all(isfinite, lb) && all(isfinite, ub))
        @warn "Skipping initialization because the variable bounds are not all finite."
    elseif initialize
        initializer = ScaledSobolIterator(lb, ub, ninit)
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
    x = copy(r.minimizer)
    prevm = minval
    iter = 1
    while iter < maxiter
        iter += 1
        ob, ineq, eq = update_surrogates!(smodel, surrogates, x)
        feasible = all(ineq .<= ctol) && all(abs.(eq) .<= ctol)
        if feasible && ob <= minval
            prevm = minval
            minval = ob
            minimizer = copy(x)
        end
        converged = feasible && abs((ob - prevm) / (prevm + ftol)) <= ftol
        converged && break
        r = optimize!(sub_workspace)
        m = ob
        x = copy(r.minimizer)
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
        smodel, optimizer.sub_alg, x0;
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
