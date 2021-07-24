
# MTS Algorithms
# MTS: Multiple Trajectory Search for Large Scale Global Optimization, 
# By [Lin-Yu Tseng and Chun Chen, 2008](https://sci2s.ugr.es/sites/default/files/files/TematicWebSites/EAMHCO/contributionsCEC08/tseng08mts.pdf)

using Random: randperm

# MTS Implementation

# Algs
struct MTSAlg <: AbstractOptimizer end
struct LS1Alg <: AbstractOptimizer end

# Options
@with_kw struct MTSOptions
    M = 100
    maxiter=200
    search_range_tol=1e-15
    n_foreground = 80
    n_local_search = 200
    n_local_search_test = 3
    n_local_search_best = 300
    BONUS1 = 10
    BONUS2 = 1 
    a_min = 0.4
    a_max = 0.5
    b_min = 0.1
    b_max = 0.3
    c_min = 0
    c_max = 1
    # Fixed parameters
    REDUCE_SEARCH_RANGE_FACTOR = 2.5
    SEARCH_RANGE_DEGERATE_FACTOR = 2
    Y1_INCR = 0.1
    Y1_DECR = 0.1
    X2_INCR = 0.2
end

@with_kw struct LS1Options
    M = 100
    maxiter=200
    search_range_tol=1e-15
    # Fixed parameters
    REDUCE_SEARCH_RANGE_FACTOR = 2.5
    SEARCH_RANGE_DEGERATE_FACTOR = 2
    Y1_INCR = 0.1
    Y1_DECR = 0.1
    X2_INCR = 0.2
    # Dummy parameters
    BONUS1 = 10
    BONUS2 = 1
end

# Workspaces
@params mutable struct MTSWorkspace <: Workspace
    model::VecModel
    x0::AbstractVector
    x::AbstractVector
    options::MTSOptions
    enable::BitVector
    improve::BitVector
    search_range::AbstractVector
    # Volitale variables
    optimal_x::AbstractVector
    optimal_ind::Int
    optimal_val::Real
end

@params mutable struct LS1Workspace <: Workspace
    model::VecModel
    x0::AbstractVector
    x::AbstractVector
    options::LS1Options
    enable::BitVector
    improve::BitVector
    search_range::AbstractVector
    # Volitale variables
    optimal_x::AbstractVector
    optimal_ind::Int
    optimal_val::Real
end

if debugging[]
    function Base.setproperty!(workspace::MTSWorkspace, name::Symbol, x)
        if name in [:optimal_x, :optimal_ind, :optimal_val]
            println("New ",  name, ": ", x)
        end
        setfield!(workspace, name, x)
    end
end

# Workspace constructors
function MTSWorkspace(model::VecModel, x0::AbstractVector, options::MTSOptions; kwargs...)
    @unpack box_min, box_max = model
    M = options.M
    # Initialize improve and serch range
    enable = trues(M)
    improve =  trues(M)
    search_range = [(box_max-box_min) ./ 2 for _ in 1:M]
    MTSWorkspace(model, x0, copy(x0), options, enable, improve, search_range, x0[1], -1, Inf)
end


function LS1Workspace(model::VecModel, x0::AbstractVector, options::LS1Options; kwargs...)
    @unpack box_min, box_max = model
    M = options.M
    # Initialize improve and serch range
    enable = trues(M)
    improve =  trues(M)
    search_range = [(box_max-box_min) ./ 2 for _ in 1:M]
    LS1Workspace(model, x0, copy(x0), options, enable, improve, search_range, x0[1], -1, Inf)
end

# Exposed workspace constructors
function Workspace(model::VecModel, optimizer::LS1Alg, x0::AbstractVector; options::LS1Options=LS1Options(), kwargs...,)
    @assert length(x0) > 0 && x0[1] isa AbstractVector
    if length(model.ineq_constraints) > 0 || length(model.eq_constraints) > 0
        @warn "LS1 does not support (in)equality constraints. Your input would be ignored. "
    end
    return LS1Workspace(model, x0, options)
end

# LS1 Workspace constructor without x0 (use method in paper to initialize)
function Workspace(model::VecModel, optimizer::LS1Alg; options::LS1Options=LS1Options(), kwargs...)
    x0 = initialize_x(model, options)
    return Workspace(model, optimizer, x0; options=options)
end

@params struct MTSResult <: AbstractResult
    minimum
    minimizer
end

# Tool functions
function initialize_x(model::VecModel, options::Union{MTSOptions, LS1Options})
    @unpack n_dim, box_min, box_max = model
    @unpack M = options
    SOA = build_SOA(M, n_dim, M)
    x0 = Array{Real,2}(undef, M, n_dim)
    for i in 1:M
        for j in 1:n_dim
            # To be confirmed: At the 5th line of "Multi Trajectory Search" algorithm in paper, I do think (u_i, l_i) shoule be (u_j, l_j)
            x0[i, j] = box_min[j] + (box_max[j]-box_min[j])*(SOA[i, j]/(M-1))
        end
    end
    [x0[i, :] for i in 1:size(x0, 1)]
end

function reduce_search_range(search_range, k, n_dim, box_min, box_max, search_range_tol, REDUCE_SEARCH_RANGE_FACTOR)
    search_range[k] ./= 2
    for i in 1:n_dim
        if search_range[k][i] < search_range_tol
            search_range[k][i] = (box_max[i] - box_min[i]) / REDUCE_SEARCH_RANGE_FACTOR
        end
    end
end

function get_buffered_clamp_and_evaluate!()
    buffer = Dict()
    function buffered_clamp_and_evaluate!(model::AbstractModel, x::AbstractVector)
        if !haskey(buffer, (model, x))      
            buffer[((model, x))] = clamp_and_evaluate!(model::AbstractModel, x::AbstractVector)
        end
        return buffer[((model, x))]
    end
    return buffered_clamp_and_evaluate!
end

# Subalgorithms
function _localsearch1(workspace, k)
    @unpack model, options = workspace
    @unpack x, improve, search_range = workspace
    @unpack BONUS1, BONUS2, search_range_tol = options
    @unpack SEARCH_RANGE_DEGERATE_FACTOR, REDUCE_SEARCH_RANGE_FACTOR = options
    @unpack box_min, box_max, n_dim = model
    grade = 0
    # search_range in paper is one-dimensional. Expand it to multidimensional. 
    if improve[k] == false
        reduce_search_range(search_range, k, n_dim, box_min, box_max, search_range_tol, REDUCE_SEARCH_RANGE_FACTOR)
    end
    improve[k] = false
    buffered_clamp_and_evaluate! = get_buffered_clamp_and_evaluate!()
    for i in 1:n_dim
        # Original value
        _xk = copy(x[k])
        xk_val = buffered_clamp_and_evaluate!(model, _xk)
        update_xki = true
        _xk[i] -= search_range[k][i]
        # Value after update
        _xk_val = buffered_clamp_and_evaluate!(model, _xk)
        # Better than current best solution
        if _xk_val < workspace.optimal_val
            grade += BONUS1
            workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
        end
        # Value stays the same
        if _xk_val == xk_val
            # Restore 
            update_xki = false
        else
            # Value degenerates
            if _xk_val > xk_val
                # Restore x_k
                _xk = copy(x[k])
                _xk[i] += search_range[k][i] / SEARCH_RANGE_DEGERATE_FACTOR
                # Current value
                _xk_val = buffered_clamp_and_evaluate!(model, _xk)
                if _xk_val < workspace.optimal_val
                    grade += BONUS1
                    workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
                end
                if _xk_val >= xk_val
                    # Restore
                    update_xki = false
                else
                    grade += BONUS2
                    improve[k] = true
                end
            else
                grade += BONUS2
                improve[k] = true
            end
        end
        if update_xki
            x[k][i] = _xk[i]
        end
    end
    return grade
end

function _localsearch2(workspace, k)
    @unpack model, options = workspace
    @unpack x, improve, search_range = workspace
    @unpack BONUS1, BONUS2, search_range_tol = options
    @unpack SEARCH_RANGE_DEGERATE_FACTOR, REDUCE_SEARCH_RANGE_FACTOR = options
    @unpack box_min, box_max, n_dim = model
    grade = 0
    # search_range in paper is one-dimensional. Expand it to multidimensional. 
    if improve[k] == false
        reduce_search_range(search_range, k, n_dim, box_min, box_max, search_range_tol, REDUCE_SEARCH_RANGE_FACTOR)
    end
    improve[k] = false
    D = zeros(Int, n_dim)
    r = zeros(Int, n_dim)
    buffered_clamp_and_evaluate! = get_buffered_clamp_and_evaluate!()
    for i in 1:n_dim
        # Original value
        xk_val = buffered_clamp_and_evaluate!(model, x[k])
        update_xk = true
        D .= rand.(Ref([-1, 1]))
        r .= rand.(Ref([0, 1, 2, 3]))
        # Value after update
        _xk = copy(x[k])
        _xk .= [r[_i] == 0 ? _xk[_i]-search_range[k][i]*D[_i] : _xk[_i] for _i in 1:length(r)]
        _xk_val = buffered_clamp_and_evaluate!(model, _xk)
        if _xk_val < workspace.optimal_val
            grade += BONUS1
            workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
        end
        # Value stays the same
        if _xk_val == xk_val
            # Restore 
            update_xk = false
        else
            # Value degenerates
            if _xk_val > xk_val
                # Restore x_k
                _xk = copy(x[k])
                _xk .= [r[_i] == 0 ? _xk[_i]+(search_range[k][i]*D[_i] / SEARCH_RANGE_DEGERATE_FACTOR) : _xk[_i] for _i in 1:length(r)]
                _xk_val = buffered_clamp_and_evaluate!(model, _xk)
                if _xk_val < workspace.optimal_val
                    grade += BONUS1
                    workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
                end
                if _xk_val >= xk_val
                    update_xk = false
                else
                    grade += BONUS2
                    improve[k] = true
                end
            else
                grade += BONUS2
                improve[k] = true
            end
        end
        if update_xk
            x[k] = _xk
        end
    end
    return grade
end

function _localsearch3(workspace, k)
    @unpack model, options = workspace
    @unpack x = workspace
    @unpack BONUS1, BONUS2 = options
    @unpack Y1_INCR, Y1_DECR, X2_INCR = options
    @unpack a_min, a_max, b_min, b_max, c_min, c_max = options
    grade = 0
    _xk = copy(x[k])
    update_xk = false
    buffered_clamp_and_evaluate! = get_buffered_clamp_and_evaluate!()
    for i in 1:model.n_dim
        # Original value
        _xk_val = buffered_clamp_and_evaluate!(model, _xk)
        update_xk = true
        # Heuristic search
        _xk_x1, _xk_y1, _xk_x2  = copy(_xk), copy(_xk), copy(_xk)
        _xk_x1[i] += Y1_INCR
        _xk_y1[i] -= Y1_DECR
        _xk_x2[i] += X2_INCR
        _xk_x1_val, _xk_y1_val, _xk_x2_val = buffered_clamp_and_evaluate!(model, _xk_x1), buffered_clamp_and_evaluate!(model, _xk_y1), buffered_clamp_and_evaluate!(model, _xk_x2)
        if _xk_x1_val < workspace.optimal_val
            grade += BONUS1
            workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
        end
        if _xk_y1_val < workspace.optimal_val
            grade += BONUS1
            workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
        end
        if _xk_x2_val < workspace.optimal_val
            grade += BONUS1
            workspace.optimal_x, workspace.optimal_ind, workspace.optimal_val = copy(_xk), k, _xk_val
        end
        D1, D2, D3 = _xk_val - _xk_x1_val, _xk_val - _xk_y1_val, _xk_val - _xk_x2_val
        grade += ((D1>0) + (D2>0) + (D3>0)) * BONUS2
        if update_xk == false
            update_xk = (D1>0) || (D2>0) || (D3>0)
        end
        a = (rand() * (a_max - a_min)) + a_min
        b = (rand() * (b_max - b_min)) + b_min
        c = (rand() * (c_max - c_max)) + c_min
        _xk[i] += a*(D1 - D2) + b*(D3 - 2*D1) + c
    end
    if update_xk
        # Update x[k]
        x[k] = _xk
    end
    return grade
end

const LOCAL_SEARCH_METHODS = [_localsearch1, _localsearch2, _localsearch3]

function mts(workspace::MTSWorkspace)
    # Multiple Trajectory Search
    @unpack options, enable = workspace
    @unpack M, n_foreground, n_local_search, n_local_search_test, n_local_search_best = options
    grade_x = [-Inf for _ in 1:M]
    for i in 1:M
        if !enable[i]
            continue
        end
        LS_testgrades = [0 for _ in 1:length(LOCAL_SEARCH_METHODS)]
        for _ in 1:n_local_search_test
            LS_testgrades = [LS_testgrades[m]+_local_search(workspace, i) for (m, _local_search) in enumerate(LOCAL_SEARCH_METHODS)]
        end
        _best_local_search = LOCAL_SEARCH_METHODS[argmax(LS_testgrades)]
        grade_x[i] = 0
        for _ in 1:n_local_search
            grade_x[i] += _best_local_search(workspace, i)
        end
    end
    for _ in 1:n_local_search_best
        _localsearch1(workspace, workspace.optimal_ind)
    end
    enable_ind = reverse(sortperm(grade_x))[begin:n_foreground]
    enable[:] .= false
    enable[enable_ind] .= true
end

function optimize!(workspace::MTSWorkspace)
    options = workspace.options
    for iter in 1:options.maxiter
        if debugging[] && iter % 50 == 0
                println("Iter ", iter, " max iter: ", options.maxiter)
        end
        mts(workspace)
    end
    MTSResult(workspace.optimal_val, workspace.optimal_x) 
end

# Build a SOA (Simultaed Orthogonal Array). Refers to section 2 of the paper posted above. 
function build_SOA(m, k, q)
    SOA = Array{Int,2}(undef, m, k)
    for c in 1:k
        q_perm = randperm(q) .- 1
        p = 0
        for r in randperm(m)
            p = (p%q)+1
            SOA[r, c] = q_perm[p]
        end
    end
    SOA
end

# mts Workspace constructor with x0
function Workspace(model::VecModel, optimizer::MTSAlg, x0::AbstractVector; options::MTSOptions=MTSOptions(), kwargs...,)
    @assert length(x0) > 0 && x0[1] isa AbstractVector
    if length(model.ineq_constraints) > 0 || length(model.eq_constraints) > 0
        @warn "MTS does not support (in)equality constraints. Your input would be ignored. "
    end
    return MTSWorkspace(model, x0, options)
end

# mts Workspace constructor without x0 (use method in paper to initialize)
function Workspace(model::VecModel, optimizer::MTSAlg; options::MTSOptions=MTSOptions(), kwargs...)
    x0 = initialize_x(model, options)
    return Workspace(model, optimizer, x0; options=options)
end

# Export localsearch1 independently
function localsearch1(workspace::Union{MTSWorkspace, LS1Workspace})
    M = workspace.options.M
    for i in 1:M
        _localsearch1(workspace, i)
    end
end

# Export LS1 independently
function optimize!(workspace::LS1Workspace)
    options = workspace.options
    for iter in 1:options.maxiter
        if debugging[] && iter % 50 == 0
                println("Iter ", iter, " max iter: ", options.maxiter)
        end
        localsearch1(workspace)
    end
    MTSResult(workspace.optimal_val, workspace.optimal_x) 
end
