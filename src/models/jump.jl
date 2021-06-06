function get_variable_info(model::JuMP.Model)
    inv_vars = OrderedDict{Int, String}()
    inv_integer = OrderedDict{Int, Bool}()
    inv_lb = OrderedDict{Int, Float64}()
    inv_ub = OrderedDict{Int, Float64}()
    inv_inits = OrderedDict{Int, Float64}()
    for (_, v) in model.obj_dict
        if v isa VariableRef
            vs = [v]
        elseif v isa AbstractArray{<:VariableRef}
            vs = v
        else
            continue
        end
        for v in vs
            ind = v.index.value
            inv_vars[ind] = JuMP.name(v)
            inv_integer[ind] = is_binary(v) || is_integer(v)
            inv_lb[ind] = has_lower_bound(v) ? lower_bound(v) : -Inf
            inv_ub[ind] = has_upper_bound(v) ? upper_bound(v) : Inf
            inv_inits[ind] = if start_value(v) !== nothing
                start_value(v)
            elseif inv_ub[ind] == Inf && inv_lb[ind] == -Inf
                0.0
            elseif inv_ub[ind] == Inf
                inv_lb[ind] + 1.0
            elseif inv_lb[ind] == -Inf
                inv_lb[ind] - 1.0
            else
                (inv_ub[ind] + inv_lb[ind]) / 2
            end
        end
    end
    sort!(inv_integer)
    sort!(inv_lb)
    sort!(inv_ub)
    sort!(inv_inits)
    lb = OrderedDict(inv_vars[k] => v for (k, v) in inv_lb)
    ub = OrderedDict(inv_vars[k] => v for (k, v) in inv_ub)
    inits = OrderedDict(inv_vars[k] => v for (k, v) in inv_inits)
    integer = OrderedDict(inv_vars[k] => v for (k, v) in inv_integer)
    return lb, ub, inits, integer
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
    for (_, v) in model.obj_dict
        if v isa ConstraintRef
            set = constraint_object(v).set
            if set isa MOI.GreaterThan
                inv_geq_constraints[JuMP.index(v).value] = v
            elseif set isa MOI.LessThan
                inv_leq_constraints[JuMP.index(v).value] = v
            elseif set isa MOI.EqualTo
                inv_eq_constraints[JuMP.index(v).value] = v
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

    for v in geq_constraints
        func = constraint_object(v).func
        @assert func isa AffExpr
        set = constraint_object(v).set
        ineqcounter += 1
        push!(bineq, -set.lower)
        for (var, val) in func.terms
            push!(Iineq, ineqcounter)
            push!(Jineq, var.index.value)
            push!(Vineq, -val)
        end
    end
    for v in leq_constraints
        func = constraint_object(v).func
        @assert func isa AffExpr
        set = constraint_object(v).set
        ineqcounter += 1
        push!(bineq, set.upper)
        for (var, val) in func.terms
            push!(Iineq, ineqcounter)
            push!(Jineq, var.index.value)
            push!(Vineq, val)
        end
    end
    for v in eq_constraints
        func = constraint_object(v).func
        @assert func isa AffExpr
        set = constraint_object(v).set
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

    return Array(Aineq), bineq, Array(Aeq), beq
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
    Aineq, bineq, Aeq, beq = get_constraint_info(model, length(lb))
    c = get_objective_info(model, length(lb))
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
