struct NLoptAlg
    alg::Symbol
end

@params struct NLoptOptions
    nt::NamedTuple
end
function NLoptOptions(;
    ftol_rel = 1e-6,
    ftol_abs = 1e-6,
    xtol_rel = 1e-6,
    xtol_abs = 1e-6,
    kwargs...,
)
    NLoptOptions(
        merge(
            NamedTuple(kwargs),
            (
                ftol_rel = ftol_rel,
                ftol_abs = ftol_abs,
                xtol_rel = xtol_rel,
                xtol_abs = xtol_abs,
            )
        )
    )
end

@params mutable struct NLoptWorkspace
    model::Model
    problem::NLopt.Opt
    x0::AbstractVector
    options::NLoptOptions
    alg::NLoptAlg
end
function NLoptWorkspace(
    model::Model, optimizer::NLoptAlg,
    x0::AbstractVector = getinit(model);
    options = NLoptOptions(), kwargs...,
)
    problem = getnlopt_problem(optimizer.alg, model, x0)
    return NLoptWorkspace(model, problem, x0, options, optimizer)
end
struct NLoptResult{X, O, P, S, A}
    minimizer::X
    minimum::O
    problem::P
    status::S
    alg::A
end

function optimize!(workspace::NLoptWorkspace)
    @unpack problem, options, x0 = workspace
    foreach(keys(options.nt)) do k
        v = options.nt[k]
        setproperty!(problem, k, v)
    end
    minf, minx, ret = NLopt.optimize(problem, x0)
    return NLoptResult(
        copy(minx), minf, problem, ret, workspace.alg,
    )
end

function Workspace(model::AbstractModel, optimizer::NLoptAlg, args...; kwargs...,)
    return NLoptWorkspace(model, optimizer, args...; kwargs...)
end

function getnlopt_problem(alg, model::Model, x0::AbstractVector)
    eq = if length(model.eq_constraints.fs) == 0
        nothing
    else
        model.eq_constraints
    end
    ineq = if length(model.ineq_constraints.fs) == 0
        nothing
    else
        model.ineq_constraints
    end
    return getnlopt_problem(
        alg,
        getobjective(model),
        ineq,
        eq,
        x0,
        getmin(model),
        getmax(model),
    )
end
function getnlopt_problem(alg, obj, ineq_constr, eq_constr, x0, xlb, xub)
    onehot(n, i) = [zeros(i-1); 1.0; zeros(n-i)]
    local lastx, objval, objgrad
    local ineqconstrval, ineqconstrjac
    if ineq_constr !== nothing
        ineqconstrval, ineqconstrjac = nothing, nothing
    end
    local eqconstrval, eqconstrjac
    if eq_constr !== nothing
        eqconstrval, eqconstrjac = nothing, nothing
    end

    function update_cached_values(x, updategrad = true)
        lastx = copy(x)
        if updategrad
            objval, objpb = Zygote.pullback(obj, lastx)
            objgrad = objpb(1.0)[1]
        else
            objval = obj(lastx)
        end
        if ineq_constr !== nothing
            if updategrad
                ineqconstrval, ineqconstrpb = Zygote.pullback(
                    ineq_constr, lastx,
                )
                ineqconstrjac = mapreduce(vcat, 1:length(ineqconstrval)) do i
                    ineqconstrpb(onehot(length(ineqconstrval), i))[1]'
                end
            else
                ineqconstrval = ineq_constr(lastx)
            end
        end
        if eq_constr !== nothing
            if length(grad) > 0
                eqconstrval, eqconstrpb = Zygote.pullback(
                    eq_constr, lastx,
                )
                eqconstrjac = mapreduce(vcat, 1:length(eqconstrval)) do i
                    eqconstrpb(onehot(length(eqconstrval), i))[1]'
                end
            else
                eqconstrval = eq_constr(lastx)
            end
        end
    end
    update_cached_values(x0)

    function nlopt_obj(x, grad)
        if x == lastx
            if length(grad) > 0
                grad .= objgrad
            end
            return objval
        else
            update_cached_values(x, length(grad) > 0)
            if length(grad) > 0
                grad .= objgrad
            end
            return objval
        end
        return val
    end
    function nlopt_ineqconstr_gen(i)
        function (x, grad)
            if x == lastx
                if length(grad) > 0
                    grad .= ineqconstrjac[i, :]
                end
                return ineqconstrval[i]
            else
                update_cached_values(x, length(grad) > 0)
                if length(grad) > 0
                    grad .= ineqconstrjac[i, :]
                end
                return ineqconstrval[i]
            end
        end
    end
    function nlopt_eqconstr_gen(i)
        function (x, grad)
            if x == lastx
                if length(grad) > 0
                    grad .= eqconstrjac[i, :]
                end
                return eqconstrval[i]
            else
                update_cached_values(x, length(grad) > 0)
                if length(grad) > 0
                    grad .= eqconstrjac[i, :]
                end
                return eqconstrval[i]
            end
        end
    end

    problem = NLopt.Opt(alg, length(x0))
    problem.lower_bounds = xlb
    problem.upper_bounds = xub

    problem.min_objective = nlopt_obj
    if ineq_constr !== nothing
        for i in 1:length(ineqconstrval)
            NLopt.inequality_constraint!(
                problem, nlopt_ineqconstr_gen(i), 0.0,
            )
        end
    end
    if eq_constr !== nothing
        for i in 1:length(eqconstrval)
            NLopt.equality_constraint!(
                problem, nlopt_eqconstr_gen(i), 0.0,
            )
        end
    end
    return problem
end
