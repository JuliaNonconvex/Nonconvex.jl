abstract type AbstractMMAApprox <: AbstractFunction end

"""
struct MMAApprox{T} <: AbstractMMAApprox
    parent::AbstractFunction # exact function
    xk::AbstractVector # current approximation point
    fk::T # f(xk)
    ∇fk::AbstractVecOrMat # ∇f(xk)
    σ::AbstractVector # move limit
    l::AbstractVector # xk .- σ
    u::AbstractVector # xk .+ σ
    ρ::T # lift - only used in MMA 2002
    p::AbstractVecOrMat # p[i, j] = σ[j]^2 max(0, ∇fk[j]) + ρ[i] σ[j] / 4
    q::AbstractVecOrMat # q[i, j] = σ[i]^2 max(0, -∇fk[j]) + ρ[i] σ[j] / 4
    r::T # r[i] = fk[i] + sum_j (p[i, j] + q[i, j]) / σ[j]
    out::T # cached output of the approximation
end

An instance of `MMAApprox` stores information about the method of moving asymptotes (MMA) approximation of the function `parent`. The notation used here is from the original [1987 MMA paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207) and the [2002 paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822). `xk` is the point around which the approximation is built. `fk` is the value of approximated function at `xk`. And `∇fk` is the gradient of the approximated function at `xk` if it is a scalar-valued function and the Jacobian if it is a vector-valued function. `l` and `u` are the lower and upper asymptotes. `σ` is the move limit such that `l = xk .- σ` and `u = xk .+ σ`. `r`, `p`, `q` and `ρ` are explained in the paper. When an instance of `MMAApprox` is called, its output will be cached in the field `out`.
"""
@params struct MMAApprox{T} <: AbstractMMAApprox
    parent::AbstractFunction # exact function
    xk::AbstractVector # current approximation point
    fk::T # f(xk)
    ∇fk::AbstractVecOrMat # ∇f(xk)
    σ::AbstractVector # move limit
    l::AbstractVector # xk .- σ
    u::AbstractVector # xk .+ σ
    ρ::T # lift - only used in MMA 2002
    p::AbstractVecOrMat # p[i, j] = σ[j]^2 max(0, ∇fk[j]) + ρ[i] σ[j] / 4
    q::AbstractVecOrMat # q[i, j] = σ[i]^2 max(0, -∇fk[j]) + ρ[i] σ[j] / 4
    r::T # r[i] = fk[i] + sum_j (p[i, j] + q[i, j]) / σ[j]
    out::T # cached output of the approximation
end

"""
    MMAApprox(f::AbstractFunction, x::AbstractVector; σ = fill(10.0, length(x)), ρ = 0.0)

Constructs an MMA approximation of the function `f` around the point `x` using a move limit `σ` and `ρ`. See [this paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822) for an explanation of the notation.
"""
function MMAApprox(f::AbstractFunction, x::AbstractVector; kwargs...)
    if getdim(f) == 1
        val, grad = value_gradient(f, x)
        return MMAApprox(f, x, val, copy(grad'); kwargs...)
    else
        val, jac = value_jacobian(f, x)
        return MMAApprox(f, x, val, jac; kwargs...)
    end
end
function MMAApprox(
    parent::AbstractFunction,
    x::AbstractVector,
    f::Real,
    ∇f::AbstractVector;
    σ = fill(10.0, length(x)),
    ρ = 0.0,
)
    l = x .- σ
    u = x .+ σ
    ρ = Ref(ρ)

    # All matrices are stored as Adjoint{<:Real, <:AbstractMatrix} for cache efficiency
    p = similar(x)
    q = similar(x)

    r = Ref(zero(f))
    out = Ref(zero(f))

    approx = MMAApprox(parent, copy(x), Ref(f), copy(∇f), σ, l, u, ρ, p, q, r, out)
    updateapprox!(approx, x, f, ∇f)
    return approx
end
function MMAApprox(
    parent::AbstractFunction,
    x::AbstractVector,
    f::Union{Real, AbstractVector},
    ∇f::AbstractMatrix;
    σ = fill(10.0, length(x)),
    ρ = f isa Real ? Ref(zero(f)) : zeros(length(f)),
)
    l = x .- σ
    u = x .+ σ

    # All matrices are stored as Adjoint{<:Real, <:AbstractMatrix} for cache efficiency
    p = similar(∇f, length(x), length(f))'
    q = similar(∇f, length(x), length(f))'

    r = f isa Real ? Ref(f) : similar(f)
    out = f isa Real ? Ref(f) : similar(f)
    fk = f isa Real ? Ref(f) : similar(f)

    approx = MMAApprox(parent, copy(x), fk, copy(∇f')', σ, l, u, ρ, p, q, r, out)
    updateapprox!(approx, x, f, ∇f)
    return approx
end

"""
    (f::MMAApprox)(x::AbstractVector)

Returns an approximation of the function `f.parent` at point `x`. See [`MMAApprox`](@ref) for an explanation.
"""
function (f::MMAApprox{<:Base.RefValue{<:Real}})(x::AbstractVector)
    @unpack p, q, l, u, r, out = f
    @assert all(l[j] < x[j] < u[j] for j in 1:length(x))
    _out = r[]
    for j in 1:length(x)
        _out += p[j] / (u[j] - x[j]) + q[j] / (x[j] - l[j])
    end
    out[] = ForwardDiff.value(_out)
    return _out
end
function (f::MMAApprox{<:AbstractVector})(x::AbstractVector)
    @unpack p, q, l, u, r, out = f
    @assert all(l[j] < x[j] < u[j] for j in 1:length(x))
    if eltype(x) <: ForwardDiff.Dual
        dual_out = map(1:getdim(f)) do i
            output = r[i]
            for j in 1:length(x)
                output += p[i, j] / (u[j] - x[j]) + q[i, j] / (x[j] - l[j])
            end
            return output
        end
        out .= ForwardDiff.value.(dual_out)
        return dual_out
    else
        for i in 1:getdim(f)
            out[i] = r[i]
            for j in 1:length(x)
                out[i] += p[i, j] / (u[j] - x[j]) + q[i, j] / (x[j] - l[j])
            end
        end
        return copy(out)
    end
end

function ChainRulesCore.rrule(f::MMAApprox{<:Base.RefValue{<:Real}}, x::AbstractVector)
    @unpack p, q, l, u, r = f
    @assert all(f.l[j] < x[j] < f.u[j] for j in 1:length(x))
    out = f(x)
    adjoint = Δ -> begin
        if p isa AbstractMatrix
            (nothing, (p' ./ (u .- x) .^ 2 .- q' ./ (x .- l) .^ 2) .* Δ)
        else
            (nothing, (p ./ (u .- x) .^ 2 .- q ./ (x .- l) .^ 2) .* Δ)
        end
    end
    f.out[] = ForwardDiff.value(out)
    return out, adjoint
end
function ChainRulesCore.rrule(f::MMAApprox{<:AbstractVector}, x::AbstractVector)
    @unpack p, q, l, u, r = f
    @assert all(f.l[j] < x[j] < f.u[j] for j in 1:length(x))
    out = f(x)
    adjoint = Δ -> begin
        (nothing, (p ./ (u' .- x') .^ 2 .- q ./ (x' .- l') .^ 2)' * Δ)
    end
    f.out .= ForwardDiff.value.(out)
    return out, adjoint
end

getmmaapprox(approx::MMAApprox) = approx

getasymptotes(approx::MMAApprox) = approx.l, approx.u

getapproxfg(approx::MMAApprox) = approx.out

getσ(approx::MMAApprox{<:Base.RefValue{<:Real}}) = approx.σ[]
getσ(approx::MMAApprox{<:AbstractVector}) = approx.σ

function setσ!(approx::MMAApprox{<:Base.RefValue{<:Real}}, σ::Real)
    approx.σ[] = σ
    return approx
end
function setσ!(approx::MMAApprox{<:AbstractVector}, σ::AbstractVector)
    approx.σ .= σ
    return approx
end

getρ(approx::MMAApprox{<:Base.RefValue{<:Real}}) = approx.ρ[]
getρ(approx::MMAApprox{<:AbstractVector}) = approx.ρ

function setρ!(approx::MMAApprox{<:Base.RefValue{<:Real}}, ρ::Real)
    approx.ρ[] = ρ
    return approx
end
function setρ!(approx::MMAApprox{<:AbstractVector}, ρ::AbstractVector)
    approx.ρ .= ρ
    return approx
end

getfk(approx::MMAApprox{<:Base.RefValue{<:Real}}) = approx.fk[]
getfk(approx::MMAApprox{<:AbstractVector}) = approx.fk

function setfk!(approx::MMAApprox{<:Base.RefValue{<:Real}}, f::Real)
    approx.fk[] = f
    return approx
end
function setfk!(approx::MMAApprox{<:AbstractVector}, f::AbstractVector)
    approx.fk .= f
    return approx
end

getxk(approx::MMAApprox) = approx.xk

function setxk!(approx::MMAApprox, x::AbstractVector)
    approx.xk .= x
    return approx
end

get∇fk(approx::MMAApprox) = approx.∇fk

function set∇fk!(approx::MMAApprox, ∇f::AbstractVecOrMat)
    approx.∇fk .= ∇f
    return approx
end

getparent(approx::MMAApprox) = approx.parent
getdim(approx::MMAApprox) = getdim(getparent(approx))

function updateapprox!(
    approx::MMAApprox{<:Base.RefValue{<:Real}},
    x::AbstractVector,
)
    val, grad = value_gradient(getparent(approx), x)
    return updateapprox!(approx, x, val, grad)
end
function updateapprox!(
    approx::MMAApprox{<:AbstractVector},
    x::AbstractVector,
)
    val, jac = value_jacobian(getparent(approx), x)
    return updateapprox!(approx, x, val, jac)
end
function updateapprox!(
    approx::MMAApprox{<:Base.RefValue{<:Real}},
    x::AbstractVector,
    f::Real,
    ∇f::Union{AbstractVector, Adjoint{<:Any, <:AbstractVector}},
)
    setxk!(approx, x)
    setfk!(approx, f)
    set∇fk!(approx, ∇f)
    return updateapprox!(approx)
end
function updateapprox!(
    approx::MMAApprox{<:AbstractVector},
    x::AbstractVector,
    f::AbstractVector,
    ∇f::AbstractMatrix,
)
    setxk!(approx, x)
    setfk!(approx, f)
    set∇fk!(approx, ∇f)
    return updateapprox!(approx)
end
function updateapprox!(approx::MMAApprox{<:Base.RefValue{<:Real}})
    @unpack xk, fk, ∇fk, l, u, σ, ρ, r, p, q, out = approx
    T = eltype(∇fk)
    l .= xk .- σ
    u .= xk .+ σ
    for j in 1:length(p)
        p[j] = σ[j] * σ[j] * max(0, ∇fk[j]) + ρ[] * σ[j] / 4
        q[j] = σ[j] * σ[j] * max(0, -∇fk[j]) + ρ[] * σ[j] / 4
    end
    r[] = fk[]
    for j in 1:length(p)
        r[] -= (p[j] + q[j]) / σ[j]
    end
    return approx
end
function updateapprox!(approx::MMAApprox{<:AbstractVector})
    @unpack xk, fk, ∇fk, l, u, σ, ρ, r, p, q, out = approx
    T = eltype(∇fk)
    l .= xk .- σ
    u .= xk .+ σ
    # All matrices are stored as Adjoint{<:Real, <:AbstractMatrix} for cache efficiency
    for i in 1:size(p, 1)
        for j in 1:size(p, 2)
            p[i, j] = σ[j] * σ[j] * max(0, ∇fk[i, j]) + ρ[i] * σ[j] / 4
            q[i, j] = σ[j] * σ[j] * max(0, -∇fk[i, j]) + ρ[i] * σ[j] / 4
        end
    end
    for i in 1:size(p, 1)
        r[i] = fk[i]
        for j in 1:size(p, 2)
            r[i] -= (p[i, j] + q[i, j]) / σ[j]
        end
    end
    return approx
end
