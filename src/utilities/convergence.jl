"""
    ConvergenceState

A struct that summarizes the convergence state of a solution. The fields in this struct are:
 - `Δx`: the infinity norm of the change in the solution `x`.
 - `Δf`: the ratio of change in the objective value `f`.
 - `kkt_residual`: the Karush-Kuhn-Tucker (KKT) residual of the solution. If [`ScaledKKTCriteria`](@ref) is used instead of [`KKTCriteria`](@ref), the `kkt_residual` will be divided by a factor.
 - `ipopt_residual`: the modified KKT residual used in the IPOPT solver.
 - `infeas`: maximum infeasibility amount. This is 0 if the solution is feasible.
 - `x_converged`: true if `Δx` is less than the `x` tolerance in [`MMAOptions`](@ref).
 - `f_converged`: true if `Δf` is less than the `f` tolerance in [`MMAOptions`](@ref).
 - `kkt_converged`: true if the `kkt_residual` is less than the KKT tolerance in [`MMAOptions`](@ref).
 - `ipopt_converged`: true if the `ipopt_residual` is less than the KKT tolerance in [`MMAOptions`](@ref).
 - `infeas_converged`: true if `infeas` is less than the infeasibility tolerance in [`MMAOptions`](@ref).
 - `f_increased`: true if the objective value of the current solution is higher than that of the previous solution.
 - `converged`: true if the solution satisfies the convergence criteria of choice. See [`ConvergenceCriteria`](@ref) for more on the different convergence criteria available.
"""
@with_kw mutable struct ConvergenceState{T}
    Δx::T = Inf
    Δf::T = Inf
    kkt_residual::T = Inf
    ipopt_residual::T = Inf
    infeas::T = Inf
    x_converged::Bool = false
    f_converged::Bool = false
    kkt_converged::Bool = false
    ipopt_converged::Bool = false
    infeas_converged::Bool = false
    f_increased::Bool = true
    converged::Bool = false
end
ConvergenceState(::Type{T}) where {T} = ConvergenceState{T}()

"""
    hasconverged(s::ConvergenceState)

Returns true if `s.converged` is true.
"""
hasconverged(s::ConvergenceState) = s.converged

"""
    Solution

A struct that stores all the information about a solution. The following are the fields of `Solution`:
 - `x`: the current primal solution
 - `prevx`: the previous primal solution
 - `λ`: the current dual solution
 - `f`: the objective value of `x`
 - `prevf`: the obejctive value of `prevx`
 - `∇f`: the gradient of the objective at `x`
 - `g`: the constraint function value at `x`
 - `∇g`: the Jacobian of the constraint functions at `x`
 - `convstate`: the convergence state of the solution
"""
@params mutable struct Solution
    prevx
    x
    λ
    prevf
    f
    ∇f
    g
    ∇g
    convstate::ConvergenceState
end

"""
    Solution(dualmodel, λ)

Construct an empty solution for the dual model `dualmodel` given a sample dual solution `λ`.
"""
function Solution(dualmodel, λ)
    prevx = copy(getxk(dualmodel))
    x = copy(prevx)
    λ = copy(λ)
    prevf = Inf
    fg = getfk(dualmodel)
    ∇fg = get∇fk(dualmodel)
    f = fg[1]
    g = fg[2:end]
    ∇f = ∇fg[1, :]
    ∇g = ∇fg[2:end, :]
    convstate = ConvergenceState()
    return Solution(prevx, x, λ, prevf, f, ∇f, g, ∇g, convstate)
end

"""
    hasconverged(s::Solution)

Returns true if `hasconverged(s.convstate)` is true.
"""
hasconverged(s::Solution) = hasconverged(s.convstate)

"""
    Tolerance

A struct specifying the different tolerances used to assess the convergence of the algorithms. The following are the fields of `Tolerance`:
 - `x`: the tolerance for `Δx` in [`ConvergenceState`](@ref). `x_converged` will be true if `Δx` is less than the `x` tolerance in `Tolerance`. This is used to assess convergence when the [`GenericCriteria`](@ref) is used as the convergence criteria. 
 - `f`: the tolerance for `Δf` in [`ConvergenceState`](@ref). `f_converged` will be true if `Δf` is less than the `f` tolerance in `Tolerance`. This is used to assess convergence when the [`GenericCriteria`](@ref) is used as the convergence criteria.
 - `kkt`: the KKT tolerance. `kkt_converged` in [`ConvergenceState`](@ref) will be true if the `kkt_residual` is less than the KKT tolerance. And `ipopt_converged` in [`ConvergenceState`](@ref) will be true if `ipopt_residual` is less than the KKT tolerance. This is used to assess convergence when the [`KKTCriteria`](@ref), the [`ScaledKKTCriteria`](@ref) or [`IpoptCriteria`](@ref) criteria is used as the convergence criteria.
 - `infeas`: the maximum infeasibility tolerance. `infeas_converged` in [`ConvergenceState`](@ref) will be true if the maximum infeasibility is less than the infeasibility tolerance. This is used to assess convergence regardless of the convergence criteria used.

For more on convergence criteria, see [`GenericCriteria`](@ref), [`KKTCriteria`](@ref), [`ScaledKKTCriteria`](@ref) and [`IpoptCriteria`](@ref).
"""
@with_kw struct Tolerance{Tx, Tf, Tkkt, Tinfeas}
    x::Tx = 0.0
    f::Tf = 1e-5
    kkt::Tkkt = 1e-5
    infeas::Tinfeas = 1e-5
end
function (tol::Tolerance{<:Function, <:Function, <:Function})(i)
    return Tolerance(tol.x(i), tol.f(i), tol.kkt(i), tol.infeas)
end

function Base.:*(t::Tolerance, m::Real)
    return Tolerance(t.x * m, t.f * m, t.kkt * m, t.infeas * m)
end

abstract type ConvergenceCriteria end

"""
    GenericCriteria

This is a generic convergence criteria that uses:
1. The maximum change in the solution, `Δx`,
2. The change percentage in the objective value, `Δf`, and
3. The maximum infeasibility `infeas`.
to assess convergence. More details are given in [`assess_convergence!`](@ref).
"""
struct GenericCriteria <: ConvergenceCriteria end

"""
    KKTCriteria

This convergence criteria uses the Karush-Kuhn-Tucker residual and maximum infeasibility to assess convergence. More details are given in [`assess_convergence!`](@ref).
"""
struct KKTCriteria <: ConvergenceCriteria end

"""
    IpoptCriteria

This convergence criteria uses a scaled version of the Karush-Kuhn-Tucker (KKT) residual and maximum infeasibility to assess convergence. This scaled KKT residual is used in the IPOPT nonlinear programming solver as explained in [this paper](http://cepac.cheme.cmu.edu/pasilectures/biegler/ipopt.pdf). More details are given in [`assess_convergence!`](@ref). 
"""
struct IpoptCriteria <: ConvergenceCriteria end

"""
    ScaledKKTCriteria

This convergence criteria uses another scaled version of the Karush-Kuhn-Tucker (KKT) residual and maximum infeasibility to assess convergence. In particular if the objective was scaled by a factor `m`, the KKT residual will be scaled down by a factor `max(m, 1/m)`. The objective scaling is explained in [`optimize!`](@ref). This scaling was found to make the convergence criteria less sensitive to scale compared to using the traditional KKT residual. More details are given in [`assess_convergence!`](@ref). 
"""
struct ScaledKKTCriteria <: ConvergenceCriteria end

"""
```
assess_convergence!(
    solution::Solution,
    model::AbstractModel,
    tol::Tolerance,
    criteria::ConvergenceCriteria,
)
```

Evaluates the convergence state `solution.convstate` given the current solution, `solution`, the tolerance, `tol`, and the convergence criteria `criteria`. `solution.convstate.converged` is then updated.

If `criteria` is an instance of `GenericCriteria`, `converged = (x_converged || f_converged) && infeas_converged`. `x_converged`, `f_converged` and `infeas_converged` are explained in [`Tolerance`](@ref). If `criteria` is an instance of `KKTCriteria` or `ScaledKKTCriteria`, `converged = kkt_converged && infeas_converged`. `kkt_converged` and `infeas_converged` are explained in [`Tolerance`](@ref). If `criteria` is an instance of `IpoptCriteria`, `converged = ipopt_converged && infeas_converged`. `ipopt_converged` and `infeas_converged` are explained in [`Tolerance`](@ref).
"""
function assess_convergence!(
    solution::Solution,
    model::AbstractModel,
    tol::Tolerance,
    criteria::ConvergenceCriteria,
)
    xtol, ftol, kkttol, infeastol = tol.x, tol.f, tol.kkt, tol.infeas
    Δx, Δf, infeas = getresiduals(solution, model, GenericCriteria())
    kkt_residual, infeas = getresiduals(solution, model, KKTCriteria())
    ipopt_residual, infeas = getresiduals(solution, model, IpoptCriteria())
    if show_residuals[]
        @show kkt_residual, ipopt_residual, kkttol
    end

    x_converged = Δx < xtol
    f_converged = Δf / (abs(solution.f) + ftol) < ftol
    if criteria isa ScaledKKTCriteria
        if debugging[]
            #@show get_objective_multiple(model)
        end
        m = get_objective_multiple(model)
        kkt_residual = kkt_residual / max(m, 1/m)
    end
    kkt_converged = kkt_residual < kkttol
    ipopt_converged = ipopt_residual < kkttol
    infeas_converged = infeas <= infeastol
    f_increased = solution.f > solution.prevf

    if criteria isa GenericCriteria
        converged = (x_converged || f_converged) && infeas_converged
    elseif criteria isa KKTCriteria || criteria isa ScaledKKTCriteria
        converged = kkt_converged && infeas_converged
    elseif criteria isa IpoptCriteria
        converged = ipopt_converged && infeas_converged
    else
        throw("Unsupported convergence criteria for MMA.")
    end
    @pack! solution.convstate = x_converged,
                                f_converged,
                                kkt_converged,
                                ipopt_converged,
                                infeas_converged,
                                Δx,
                                Δf,
                                kkt_residual,
                                ipopt_residual,
                                infeas,
                                f_increased,
                                converged
    return solution
end

function getresiduals(solution::Solution, model::AbstractModel, ::GenericCriteria)
    @unpack prevx, x, prevf, f, g = solution
    Δx = maximum(abs(x[j] - prevx[j]) for j in 1:length(x))
    Δf = abs(f - prevf)
    infeas = length(g) == 0 ? zero(eltype(g)) : max(0, maximum(g))
    return Δx, Δf, infeas
end
function getresiduals(solution::Solution, model::AbstractModel, ::KKTCriteria)
    @unpack ∇f, g, ∇g, λ, x = solution
    xmin, xmax = getmin(model), getmax(model)
    T = eltype(x)
    res = maximum(1:length(x)) do j
        @views temp = ∇f[j] + dot(∇g[:,j], λ)
        if xmin[j] >= x[j]
            return abs(min(0, temp))
        elseif x[j] >= xmax[j]
            return max(0, temp)
        else
            return abs(temp)
        end
    end
    if debugging[]
        @show λ, g
    end
    res = length(g) == 0 ? res : max(res, maximum(abs.(λ .* g)))
    if debugging[]
        @show maximum(abs, g)
        @show maximum(abs, λ)
        @show maximum(x)
    end
    infeas = length(g) == 0 ? zero(eltype(g)) : max(maximum(g), 0)
    return res, infeas
end
function getresiduals(solution::Solution, model::AbstractModel, ::IpoptCriteria)
    @unpack ∇f, g, ∇g, λ, x = solution
    xmin, xmax = getmin(model), getmax(model)
    T = eltype(x)
    n, m, s = length(x), length(λ), zero(T)
    res = maximum(1:n) do j
        @views temp = ∇f[j] + dot(∇g[:,j], λ)
        if xmin[j] >= x[j]
            dj = temp
            s += max(dj, 0)
            return abs(min(0, dj))
        elseif x[j] >= xmax[j]
            yj = -temp
            s += max(yj, 0)
            return abs(min(0, yj))
        else
            return abs(temp)
        end
    end
    sd = max(100, (sum(abs, λ) + s) / (n + m)) / 100
    res = length(g) == 0 ? res : max(res, maximum(abs.(λ .* g)))
    res = res / sd
    infeas = length(g) == 0 ? zero(eltype(g)) : max(maximum(g), 0)
    if debugging[]
        println("Agg infeas = ", infeas)
    end
    return res, infeas
end
