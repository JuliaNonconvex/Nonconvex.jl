# Semidefinite programming

# Decompress
function lowertriangind(mat::Matrix)
    indices = [i for i in CartesianIndices(mat) if i[1]>i[2]]
    return LinearIndices(mat)[indices]
end

function rearrange_x(x_L::AbstractVector, x_D::AbstractVector)
    mat_dim = length(x_D)
    L = zeros(mat_dim, mat_dim)
    L[lowertriangind(L)] .= x_L
    D = diagm(x_D)
    return (L, D)
end

"""
    decompress_symmetric

    For example: 
    a 3*3 positive semidefinite matrix: 
    `[  a       b       d; 
        b       c       e; 
        d       e       f       ]` 
    represents by: 
    `[  x_D[1]  x_L[3]  x_L[2]; 
        x_L[1]  x_D[2]  x_L[1]; 
        x_L[2]  x_L[3]  x_D[3]  ]`

- `x_L::AbstractArray`: representing lower triangular part of a `n*n` matrix, length should be `(n^2-n)÷2`
- `x_D::AbstractArray`: representing diagonal part of a `n*n` matrix, length should be `n`
"""
function decompress_symmetric(L::Matrix, D::Matrix)
    L + D + L'
end

function decompress_symmetric(x_L::AbstractArray, x_D::AbstractArray)
    L, D = rearrange_x(x_L, x_D)
    return decompress_symmetric(L, D)
end

function ChainRulesCore.rrule(::typeof(rearrange_x), x_L::AbstractVector, x_D::AbstractVector)
    function pullback((ΔL, ΔD))
        Δx_L = ΔL[lowertriangind(ΔL)]
        Δx_D = diag(ΔD)
        NoTangent(), Δx_L, Δx_D 
    end
    return rearrange_x(x_L, x_D), pullback
end


"""
    SDPBarrierAlg: a meta-algorithm used for Semidefinite barrier programming

- `sub_alg`: algorithm to solve underlying optimization objective, which should be in consistent with `sub_options` of `SDPBarrierOptions`
"""
struct SDPBarrierAlg{Alg <: AbstractOptimizer} <: AbstractOptimizer
    sub_alg::Alg
end
function SDPBarrierAlg(;sub_alg)
    return SDPBarrierAlg(sub_alg)
end


"""
    SDPBarrierOptions: options for sdp programming

- `c_init`: initial multiplier `c` for barrier, could be a real number or vector. Vector input is for the case that objective contains multiple matrix and each element of vector for corresponding one of them
- `c_decr`: decreasing rate that multiplies `c` in every iteration. Could be either a real number of a vector, same as `c_init`
- `n_iter`: how many iterations in total for barrier method
- `sub_options`: options to solve underlying optimization target
- `keep_all`: if set to `True`, `SDPBarrierResult` keeps all results in each iteration but not only optima
- `sub_alg`: algorithm to solve underlying optimization objective, which should be in consistent with `sub_options` in `SDPBarrierOptions`
"""
# Options
@params mutable struct SDPBarrierOptions
    # Dimension of objective matrix
    # Hyperparameters 
    # Initial value of `c` in barrier method: 
    # `Real` for using same value for all `sd_constraints`, `AbstractArray` for assign them respectively
    c_init::Union{Real, AbstractArray}
    # Decrease rate of `c` for every epoch, same as above
    c_decr::Union{Real, AbstractArray}
    n_iter::Int
    # sub_option to solve (in)equality constraints
    sub_options
    # Keep all results or not
    keep_all::Bool
end
function SDPBarrierOptions(c_init, c_decr, n_iter; sub_options, keep_all=true)
    @assert 0 < c_decr < 1 "c_decr should be between 0 and 1. "
    @assert c_init > 0 "c_init shoule be larger than 0. "
    SDPBarrierOptions(c_init, c_decr, n_iter, sub_options, keep_all)
end
function SDPBarrierOptions(;sub_options, c_init=1.0, c_decr=0.1, n_iter=10, keep_all=true)
    SDPBarrierOptions(c_init, c_decr, n_iter, sub_options=sub_options, keep_all=keep_all)
end

# Result
@params struct SDPBarrierResult
    minimum
    minimizer
    results
    optimal_ind
end

@params struct SDPBarrierWorkspace <: Workspace
    model::VecModel
    x0::AbstractVector
    options::SDPBarrierOptions
    sub_alg::AbstractOptimizer
end

function Workspace(model::VecModel, optimizer::SDPBarrierAlg, x0, args...; options, kwargs...,)
    @unpack c_init, c_decr = options
    for c in model.sd_constraints.fs
        @assert isposdef(c(x0)) "Initial matrix should be positive definite. "
    end
    if c_init isa AbstractArray
        @assert length(model.sd_constraints.fs) == length(c_init) "c_init should be same length with number of `sd_constraints` when using array. "
    end
    if c_decr isa AbstractArray
        @assert length(model.sd_constraints.fs) == length(c_decr) "c_decr should be same length with number of `sd_constraints` when using array. "
    end
    if c_init isa AbstractArray && c_decr isa AbstractArray
        @assert length(c_init) == length(c_decr) "c_decr should be same length with c_init. "
    end
    return SDPBarrierWorkspace(model, copy(x0), options, optimizer.sub_alg)
end

function sd_objective(objective0, sd_constraints, c::AbstractArray)
    function _objective(args)
        target = objective0(args)
        barrier = sum(c .* -logdet.(map(f -> f(args), sd_constraints.fs)))
        return target + barrier
    end
    return _objective
end

function to_barrier(model, c::AbstractArray)
    sd_constraints, objective0 = model.sd_constraints, model.objective
    _model = set_objective(model, sd_objective(objective0, sd_constraints, c))
    return _model
end

function optimize!(workspace::SDPBarrierWorkspace)
    @unpack model, x0, options, sub_alg = workspace
    @unpack c_init, c_decr, n_iter, sub_options, keep_all = options
    objective0 = model.objective
    x = copy(x0)
    c = c_init isa Real ? ([c_init for _ in 1:length(model.sd_constraints.fs)]) : c_init
    results = []
    for _ in 1:n_iter
        model_i = to_barrier(model, c)
        result_i = optimize(model_i, sub_alg, x, options = sub_options)
        minimizer_i = result_i.minimizer
        push!(results, (objective0(minimizer_i), minimizer_i))
        c = c .* c_decr
        x = copy(minimizer_i)
    end
    optimal_ind = argmin(first.(results))
    minimum, minimizer = results[optimal_ind]
    if keep_all
        return SDPBarrierResult(minimum, minimizer, results, optimal_ind)
    else
        return SDPBarrierResult(minimum, minimizer, nothing, nothing)
    end
end
