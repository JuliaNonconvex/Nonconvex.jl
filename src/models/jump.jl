function get_variable_info(model::JuMP.Model)
    varinds = OrderedDict()
    integer = OrderedDict()
    lb = OrderedDict()
    ub = OrderedDict()
    inits = OrderedDict()
    for (k, vs) in model.obj_dict
        if vs isa VariableRef
            v = vs
            ind = v.index.value
            varinds[k] = ind
            integer[k] = is_binary(v) || is_integer(v)
            lb[k] = has_lower_bound(v) ? lower_bound(v) : -Inf
            ub[k] = has_upper_bound(v) ? upper_bound(v) : Inf
            inits[k] = if start_value(v) !== nothing
                start_value(v)
            elseif ub[k] == Inf && lb[k] == -Inf
                0.0
            elseif ub[k] == Inf
                lb[k] + 1.0
            elseif lb[k] == -Inf
                lb[k] - 1.0
            else
                (ub[k] + lb[k]) / 2
            end
        elseif vs isa AbstractArray{<:VariableRef}
            inds = map(v -> v.index.value, vs)
            varinds[k] = collect(keys(inds))
            integer[k] = map(vs) do v
                is_binary(v) || is_integer(v)
            end
            lb[k] = map(vs) do v
                has_lower_bound(v) ? lower_bound(v) : -Inf
            end
            ub[k] = map(vs) do v
                has_upper_bound(v) ? upper_bound(v) : Inf
            end
            inits[k] = map(vs) do v
                ind = varinds[k][v.index.value]
                if start_value(v) !== nothing
                    start_value(v)
                elseif ub[k][ind] == Inf && lb[k][ind] == -Inf
                    0.0
                elseif ub[k][ind] == Inf
                    lb[k][ind] + 1.0
                elseif lb[k][ind] == -Inf
                    lb[k][ind] - 1.0
                else
                    (ub[k][ind] + lb[k][ind]) / 2
                end
            end
        else
            continue
        end
    end
    sort!(integer)
    sort!(lb)
    sort!(ub)
    sort!(inits)
    _lb = OrderedDict(k => v for (k, v) in lb)
    _ub = OrderedDict(k => v for (k, v) in ub)
    _inits = OrderedDict(k => v for (k, v) in inits)
    _integer = OrderedDict(k => v for (k, v) in integer)
    return _lb, _ub, _inits, _integer
end

function get_constraint_info(model::JuMP.Model, nvars)
    Iineq = Int[]
    Jineq = Int[]
    Vineq = Float64[]
    bineq = Float64[]

    Ieq = Int[]
    Jeq = Int[]
    Veq = Float64[]
    beq = Float64[]

    ineqcounter = 0
    eqcounter = 0
    inv_geq_constraints = OrderedDict()
    inv_leq_constraints = OrderedDict()
    inv_eq_constraints = OrderedDict()
    for (_, c) in model.obj_dict
        if c isa ConstraintRef
            set = constraint_object(c).set
            if set isa MOI.GreaterThan
                inv_geq_constraints[JuMP.index(c).value] = c
            elseif set isa MOI.LessThan
                inv_leq_constraints[JuMP.index(c).value] = c
            elseif set isa MOI.EqualTo
                inv_eq_constraints[JuMP.index(c).value] = c
            else
                throw("Set type not supported.")
            end
        end
    end
    sort!(inv_geq_constraints)
    sort!(inv_leq_constraints)
    sort!(inv_eq_constraints)
    geq_constraints = values(inv_geq_constraints)
    leq_constraints = values(inv_leq_constraints)
    eq_constraints = values(inv_eq_constraints)

    for c in geq_constraints
        func = constraint_object(c).func
        @assert func isa AffExpr
        set = constraint_object(c).set
        ineqcounter += 1
        push!(bineq, -set.lower)
        for (var, val) in func.terms
            push!(Iineq, ineqcounter)
            push!(Jineq, var.index.value)
            push!(Vineq, -val)
        end
    end
    for c in leq_constraints
        func = constraint_object(c).func
        @assert func isa AffExpr
        set = constraint_object(c).set
        ineqcounter += 1
        push!(bineq, set.upper)
        for (var, val) in func.terms
            push!(Iineq, ineqcounter)
            push!(Jineq, var.index.value)
            push!(Vineq, val)
        end
    end
    for c in eq_constraints
        func = constraint_object(c).func
        @assert func isa AffExpr
        set = constraint_object(c).set
        eqcounter += 1
        push!(beq, set.value)
        for (var, val) in func.terms
            push!(Ieq, eqcounter)
            push!(Jeq, var.index.value)
            push!(Veq, val)
        end
    end
    Aineq = sparse(Iineq, Jineq, Vineq, ineqcounter, nvars)
    Aeq = sparse(Ieq, Jeq, Veq, eqcounter, nvars)

    return Aineq, bineq, Aeq, beq
end

function get_objective_info(model::JuMP.Model, nvars)
    obj = objective_function(model)
    @assert obj isa AffExpr
    sense = objective_sense(model)
    I = Int[]
    V = Float64[]
    for (var, val) in obj.terms
        push!(I, var.index.value)
        if sense == MathOptInterface.MIN_SENSE
            push!(V, val)
        elseif sense == MathOptInterface.MAX_SENSE
            push!(V, -val)
        end
    end
    c = sparsevec(I, V, nvars)
    return c
end

function DictModel(model::JuMP.Model)
    lb, ub, inits, integer = get_variable_info(model)
    nvars = length(flatten(lb)[1])
    Aineq, bineq, Aeq, beq = get_constraint_info(model, nvars)
    c = get_objective_info(model, nvars)
    ks = collect(keys(inits))
    dict_model = DictModel()
    for k in ks
        addvar!(dict_model, k, lb[k], ub[k], init = inits[k], integer = integer[k])
    end
    if length(c) > 0
        set_objective!(dict_model, x -> dot(c, flatten(x, ks)[1]))
    end
    if length(bineq) > 0
        add_ineq_constraint!(dict_model, x -> Aineq * flatten(x, ks)[1] - bineq)
    end
    if length(beq) > 0
        add_eq_constraint!(dict_model, x -> Aeq * flatten(x, ks)[1] - beq)
    end
    return dict_model
end
