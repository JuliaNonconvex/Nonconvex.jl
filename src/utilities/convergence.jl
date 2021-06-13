
"""
hasconverged(s::ConvergenceState)

Returns true if `s.converged` is true.
"""
hasconverged(s::ConvergenceState) = s.converged

"""
hasconverged(s::Solution)

Returns true if `hasconverged(s.convstate)` is true.
"""
hasconverged(s::Solution) = hasconverged(s.convstate)


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
    xtol, fabstol, freltol, kkttol, infeastol = tol.x, tol.fabs, tol.frel, tol.kkt, tol.infeas
    Δx, Δf, infeas = getresiduals(solution, model, GenericCriteria())
    relΔf = Δf / (abs(solution.f) + freltol)
    kkt_residual, infeas = getresiduals(solution, model, KKTCriteria())
    ipopt_residual, infeas = getresiduals(solution, model, IpoptCriteria())
    if show_residuals[]
        @show kkt_residual, ipopt_residual, kkttol
    end

    x_converged = Δx < xtol
    fabs_converged = Δf < fabstol
    frel_converged = relΔf < freltol
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
        converged = (x_converged || fabs_converged || frel_converged) && infeas_converged
    elseif criteria isa KKTCriteria || criteria isa ScaledKKTCriteria
        converged = kkt_converged && infeas_converged
    elseif criteria isa IpoptCriteria
        converged = ipopt_converged && infeas_converged
    else
        throw("Unsupported convergence criteria for MMA.")
    end
    @pack! solution.convstate = x_converged,
                                fabs_converged,
                                frel_converged,
                                kkt_converged,
                                ipopt_converged,
                                infeas_converged,
                                Δx,
                                Δf,
                                relΔf,
                                kkt_residual,
                                ipopt_residual,
                                infeas,
                                f_increased,
                                converged
    return solution
end

function getresiduals(solution::Solution, ::AbstractModel, ::GenericCriteria)
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
