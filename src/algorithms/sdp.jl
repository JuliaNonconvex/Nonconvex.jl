# Semidefinite programming

# Decompress
function rearrange_x(x_L::AbstractVector, x_D::AbstractVector)
    mat_dim = length(x_D)
    L = zeros(mat_dim, mat_dim)
    L[lowertriangind(L)] .= x_L
    D = diagm(x_D)
    return (L, D)
end

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


# Optimizer (algorithm)
struct SDPBarrierAlg <: AbstractOptimizer
    sub_alg::AbstractOptimizer
end
function SDPBarrierAlg(;sub_alg)
    return SDPBarrierAlg(sub_alg)
end

# Options
mutable struct SDPBarrierOptions
    # Dimension of objective matrix
    # Hyperparameters 
    c_init::Real
    c_decr::Real
    n_iter::Int
    # Dimension of input matrix
    mat_dim::Int
    # sub_option to solve (in)equality constraints
    sub_options
    # Semidefinite function, returns target
    # criteria to keep semidefiniteness
    sd_criteria::String
    # Keep all results or not
    keep_all::Bool
end
function SDPBarrierOptions(c_init, c_decr, n_iter, mat_dim; sub_options, sd_criteria="neg_determinant", keep_all=true)
    @assert 0 < c_decr < 1 "c_decr should be between 0 and 1. "
    @assert c_init > 0 "c_init shoule be larger than 0. "
    SDPBarrierOptions(c_init, c_decr, n_iter, mat_dim, sub_options, sd_criteria, keep_all)
end
function SDPBarrierOptions(;c_init, c_decr, n_iter, mat_dim, sub_options, sd_criteria="neg_determinant", keep_all=true)
    SDPBarrierOptions(c_init, c_decr, n_iter, mat_dim, sub_options=sub_options, sd_criteria=sd_criteria, keep_all=keep_all)
end

# Result
@params struct SDPResult
    minimum
    minimizer
    results
    optimal_ind
end

@params struct SDPWorkspace <: Workspace
    model::VecModel
    x0::AbstractVector
    options::SDPBarrierOptions
    sub_alg::AbstractOptimizer
end

function Workspace(model::VecModel, optimizer::SDPBarrierAlg, x0, args...; options, kwargs...,)
    A = model.sd_function(x0)
    @assert all(eigvals(A) .≥ 0) "x0 matrix should be positive semidefinite. "
    return SDPWorkspace(model, copy(x0), options, optimizer.sub_alg)
end

_sd_criterias = Dict(
    "neg_determinant" => (mat) -> -logdet(mat)
)

function sd_objective(objective0, sd_function, sd_criteria, c)
    function _objective(args)
        target = objective0(args)
        barrier = c * _sd_criterias[sd_criteria](sd_function(args))
        return target + barrier
    end
    return _objective
end

function to_barrier(model::VecModel, c::Real, sd_criteria::String)
    @unpack objective, sd_function = model
    _model = deepcopy(model)
    set_objective!(_model, sd_objective(objective, sd_function, sd_criteria, c))
    return _model
end

function optimize!(workspace::SDPWorkspace)
    @unpack model, x0, options, sub_alg = workspace
    @unpack c_init, c_decr, n_iter, mat_dim, sub_options, sd_criteria, keep_all = options
    c = c_init
    results = []
    for _ in 1:n_iter
        model_i = to_barrier(model, c, sd_criteria)
        result_i = optimize(model_i, sub_alg, copy(x0), options = sub_options)
        push!(results, (result_i.minimum, result_i.minimizer))
        c *= c_decr
    end
    optimal_ind = argmin(first.(results))
    minimum, minimizer = results[optimal_ind]
    if keep_all
        return SDPResult(minimum, minimizer, results, optimal_ind)
    else
        return SDPResult(minimum, minimizer, nothing, nothing)
    end
end
