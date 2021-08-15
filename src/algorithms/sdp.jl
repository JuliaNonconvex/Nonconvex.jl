# Semidefinite programming

# Types
MatTypeA = Union{Matrix}
TypeA = Union{Real, Matrix}

# Util functions
function dof_to_dim(dof::Int)
    @assert is_triangular(dof) "dof must be a triangular number. "
    return trunc(Int, (sqrt(8*dof+1)-1)÷2)
end

function dof_to_length(dof::Int)
    dim = dof_to_dim(dof)
    return dim^2
end

dim_to_dof(dim::Int) = (dim^2+dim)÷2

function apply_sd_constraint(A0::TypeA, Ai::AbstractVector, x::AbstractVector)
    return reduce(.+, Ai.*x).+A0
end

function rearrange_x(x::AbstractVector)
    d = trunc(Int, (sqrt(8*length(x)+1)-1)÷2)
    size_L = (d^2-d)÷2
    size_D = d
    D = diagm(x[size_L+1:size_L+size_D])
    L = diagm(ones(d))
    L[lowertriangind(L)] .= x[begin:size_L]
    return (L, D)
end

function ChainRulesCore.rrule(::typeof(rearrange_x), x::AbstractVector)
    function pullback(Δ)
        ΔL, ΔD = Δ
        Δx = [ΔL[lowertriangind(ΔL)]..., diag(ΔD)...]
        NoTangent(), Δx
    end
    return rearrange_x(x), pullback
end


function inv_ldl(L::MatTypeA, D::MatTypeA; flatten=true)
    M = L*D*L'
    return flatten ? vec(M) : M
end

function inv_ldl(x::AbstractVector; flatten=true)
    L, D = rearrange_x(x)
    inv_ldl(L, D; flatten=flatten)
end

# Optimizer (algorithm)
struct SDPAlg <: AbstractOptimizer
    sub_alg::AbstractOptimizer
end
function SDPAlg(;sub_alg)
    return SDPAlg(sub_alg)
end

# Options
@params mutable struct SDPOptions
    # Dimension of objective matrix 
    mat_dim::Int
    sub_options
    sd_criteria
    # Semidefinite constraint
    # Constant constraint
    A0::TypeA
    # The first (obj_dim - dof + mat_length) elements are for matrix after inverse LDL
    # And elements come after are for `x[dof+1:end]`, if there are x remaining (`|Ai|<obj_dim - dof + mat_length`), assign remaining elements with 0
    Ai::AbstractVector
end
function SDPOptions(mat_dim; sub_options, sd_criteria="neg_determinant", A0=0, Ai=[])
    SDPOptions(mat_dim, sub_options, sd_criteria, A0, Ai)
end

function set_A0(options::SDPOptions, A0::TypeA)
    options.A0 = A0
end

function add_Ai(options::SDPOptions, A::TypeA)
    # first_mat = findfirst(e -> e isa Matrix, options._A)
    @unpack A0, Ai = options
    first_mat = A0 isa Matrix ? A0 : length(Ai) == 0 ? nothing : Ai[findfirst(e -> e isa MatTypeA, Ai)]
    @assert first_mat isa Nothing || size(first_mat) == size(A) "Semidefinite constraints should all be same size. "
    push!(options.Ai, A)
end

get_mat_dim(options::SDPOptions) = options.mat_dim
get_mat_length(options::SDPOptions) = get_mat_dim(options)^2

function get_dof(options::SDPOptions) 
    d = get_mat_dim(options)
    return dim_to_dof(d)
end

function get_length_after_transformation(obj_dim, options::SDPOptions)
    dof, mat_length = get_dof(options), get_mat_length(options)
    obj_dim - dof + mat_length
end

# Workspace
const _SD_CRITERIAS = Dict(
    "neg_determinant" => (mat) -> -det(mat)
)

# Result
@params struct SDPResult
    sub_result::AbstractResult
    minimizer_mat::Matrix
    minimum
    minimizer
end
function SDPResult(sub_result::AbstractResult, options::SDPOptions)
    @unpack minimum, minimizer = sub_result
    mat_length = get_mat_length(options)
    minimizer_mat = inv_ldl(minimizer[begin:mat_length], flatten=false)
    return SDPResult(sub_result, minimizer_mat, minimum, minimizer)
end

# Transform x from LD to matrix, then concatenate remaining elements
inv_ldl_x(x::AbstractVector, dof::Int) = vcat(inv_ldl(x[begin:dof], flatten=true), x[dof+1:end])

@params struct SDPWorkspace <: Workspace
    sub_workspace::Workspace
    options::SDPOptions
    x0::AbstractVector
end

function Workspace(model::VecModel, alg::SDPAlg, x0::AbstractVector, args...; options::SDPOptions, kwargs...)
    # Create sub-workspace
    sub_workspace = Workspace(model, alg.sub_alg, x0, args...; options=options.sub_options, kwargs...)
    return SDPWorkspace(sub_workspace, options, copy(x0))
end

function optimize!(workspace::SDPWorkspace)
    sub_result = optimize!(workspace.sub_workspace)
    return SDPResult(sub_result, workspace.options)
end

function pre_tovecmodel!(model::Model, optimizer::SDPAlg, x0, args...; options, kwargs...)
    mat_f = model.objective.f
    obj_f = mat_f.f
    # Check if input is wrapped
    @assert mat_f isa MatrixFunctionWrapper "Please wrap function using `MatrixFunctionWrapper` when performing SDP. "
    @unpack sd_criteria, A0, Ai, sub_options = options
    mat_dim, obj_dim = mat_f.mat_dim, obj_f.dim
    # Construct matrix
    dof = get_dof(options)
    # Check if x0 has correct length
    @assert obj_dim == length(x0) "Objective dim should be equal to length of x0. " 
    # Check if input dim >= dof
    @assert obj_dim >= dof "Input dim should be at least equal to dof (degrees of freedom) of matrix, dof=(mat_dim^2+mat_dim)÷2. "
    # Check if length of transformed input is at least equals to length of Ai
    @assert get_length_after_transformation(obj_dim, options) >= length(Ai) "Length of transformed input is at least equals to length of Ai. "
    # Transform model objective
    set_objective!(model, x->model.objective(inv_ldl_x(x, dof)))
    # Transfer semidefinite constraints to ineq constraints
    add_ineq_constraint!(model, (x -> _SD_CRITERIAS[sd_criteria](apply_sd_constraint(A0, Ai, inv_ldl(x[begin:dof], flatten=true)))))
end
