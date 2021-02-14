"""
```
struct XMMAApprox <: AbstractMMAApprox
    mma::MMAApprox
    approxfg::AbstractVector
    a::AbstractVector
    c::AbstractVector
    d::AbstractVector
    y::AbstractVector
    z::Base.RefValue{<:Real}
    out::AbstractVector
end
```

This is similar to [`MMAApprox`](@ref) but it uses the extended formulation given in the [2002 MMA paper](https://epubs.siam.org/doi/abs/10.1137/S1052623499362822). `y` and `z` are additional variables added and `a`, `c` and `d` are coefficients, all of which are explained in the paper. `mma` is the underlying MMA approximation of the objective and constraints. Calling `mma` returns a vector whose first element is the objective approximation and the rest are the constraint violation approximations.

If the original model is:
 - Variables: x
 - Objective: f(x)
 - Constraint i: g(x) <= 0
 - Bounds: x ∈ X
the transformed model is:
 - Variables: x, y, z
 - Objective: f(x) + a[1] * z + Σ_i c[i] * y[i] + 1/2 * d[i] * y[i]^2
 - Constraint i: g_i(x) - a[i + 1] * z - y[i] <= 0
 - Bounds: x ∈ X, y >= 0, z >= 0
"""
@params struct XMMAApprox <: AbstractMMAApprox
    mma::MMAApprox
    approxfg::AbstractVector
    a::AbstractVector
    c::AbstractVector
    d::AbstractVector
    y::AbstractVector
    z::Base.RefValue{<:Real}
    out::AbstractVector
end

"""
```
XMMAApprox(
    mma::MMAApprox,
    a = fill(1.0, size(mma.p, 1)),
    c = fill(1.1, size(mma.p, 1) - 1),
    d = fill(1.0, size(mma.p, 1) - 1),
    kwargs...,
)
```

Constructs an instance of the extended MMA approximation `XMMAApprox` using the MMA approximation `mma` and coefficients `a`, `b` and `c`. See [`XMMAApprox`](@ref) for an explanation.
"""
function XMMAApprox(
    mma::MMAApprox,
    a = fill(1.0, size(mma.p, 1)),
    c = fill(1.1, size(mma.p, 1) - 1),
    d = fill(1.0, size(mma.p, 1) - 1),
    kwargs...,
)
    x0 = getxk(mma)
    T = eltype(x0)
    approxfg = mma(x0)
    z0 = sqrt(eps(T))
    y0 = max.(approxfg[2:end] .- a[2:end] .* z0, z0)
    out = similar(mma.out)
    return XMMAApprox(mma, approxfg, a, c, d, y0, Ref(z0), out)
end

"""
    (approx::XMMAApprox)(xyz::AbstractVector)

Returns the values of the approximate objective and constraints. `xyz` is `x`, `y` and `z` concatenated in a single vector.
"""
function (approx::XMMAApprox)(xyz::AbstractVector)
    @unpack mma, a, c, d = approx
    T = eltype(xyz)
    z = xyz[end]
    y = xyz[end - length(c) : end - 1]
    x = xyz[begin : end - length(c) - 1]
    approxfg = mma(x)
    ε = 1e-4
    out = approxfg .+ [a[1] * z + ε * z^2 + sum(c .* y .+ d .* y.^2 ./ 2); .-a[2:end] .* z .- y]
    saveout!(approx, out)
    return out
end

"""
    saveout!(approx::XMMAApprox, out)

Saves `out` in the cache `approx.out`.
"""
function saveout!(approx::XMMAApprox, out)
    approx.out .= out
    return approx
end
ChainRulesCore.@non_differentiable saveout!(approx::XMMAApprox, out)

getmmaapprox(approx::XMMAApprox) = approx.mma

getasymptotes(approx::XMMAApprox) = getasymptotes(approx.mma)

getapproxfg(approx::XMMAApprox) = approx.mma.out

getσ(approx::XMMAApprox) = getσ(approx.mma)
function setσ!(approx::XMMAApprox, σ::AbstractVector)
    setσ!(approx.mma, σ)
    return approx
end

getρ(approx::XMMAApprox) = getρ(approx.mma)
function setρ!(approx::XMMAApprox, ρ::AbstractVector)
    setρ!(approx.mma, ρ)
    return approx
end

getfk(approx::XMMAApprox) = getfk(approx.mma)
function setfk!(approx::XMMAApprox, f::AbstractVector)
    setfk!(approx.mma, f)
    return approx
end

getxk(approx::XMMAApprox) = getxk(approx.mma)
function setxk!(approx::XMMAApprox, x::AbstractVector)
    setxk!(approx.mma, x)
    return approx
end

get∇fk(approx::XMMAApprox) = get∇fk(approx.mma)
function set∇fk!(approx::XMMAApprox, ∇f::AbstractVecOrMat)
    set∇fk!(approx.mma, ∇f)
    return approx
end

getparent(approx::XMMAApprox) = getparent(approx.mma)
getdim(approx::XMMAApprox) = getdim(approx.mma)

function updateapprox!(
    approx::XMMAApprox,
    x::AbstractVector,
)
    updateapprox!(approx.mma, x)
    return approx
end
function updateapprox!(
    approx::XMMAApprox,
    x::AbstractVector,
    f::AbstractVector,
    ∇f::AbstractMatrix,
)
    updateapprox!(approx.mma, x, f, ∇f)
    return approx
end
function updateapprox!(approx::XMMAApprox)
    updateapprox!(approx.mma)
    return approx
end
