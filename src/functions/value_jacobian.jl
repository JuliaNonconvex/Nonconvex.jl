# Simultaneous function and Jacobian evaluation 

"""
    value_jacobian(f::AbstractFunction, x::AbstractVector)

Returns the value and Jacobian of the function `f` at the point `x`. The automatic differentiation package, Zygote, is used by default. To define a custom adjoint rule for a function to be used in a constraint or objective, use [`ChainRulesCore.jl`](https://github.com/JuliaDiff/ChainRulesCore.jl). The Jacobian is returned as an instance of `Adjoint` for cache efficiency in the optimizer.
"""
function value_jacobian(f::AbstractFunction, x::AbstractVector)
    out, pullback = Zygote.pullback(f, x)
    @assert length(out) == getdim(f)
    if out isa Number
        grad = pullback(1.0)[1]
        return out, grad'
    else
        dim = getdim(f)
        jact = mapreduce(hcat, 1:dim; init = zeros(length(x), 0)) do i
            pullback(I(dim)[:,i])[1]
        end
        return out, jact'
    end
end
function value_jacobian(c::VectorOfFunctions, x::AbstractVector)
    vals_jacs = map(c.fs) do f
        value_jacobian(f, x)
    end
    vals = mapreduce(vcat, vals_jacs; init = zeros(0)) do val_jac
        val_jac[1]
    end
    jact = mapreduce(hcat, vals_jacs; init = zeros(length(x), 0)) do val_jac
        val_jac[2]'
    end
    return vals, jact'
end

"""
    value_jacobian_transpose(f::AbstractFunction, x::AbstractVector)

Returns the value and transpose of the Jacobian of `f` at the point `x`. This calls [`value_jacobian`](@ref) and transposes the Jacobian.
"""
function value_jacobian_transpose(f::AbstractFunction, x::AbstractVector)
    val, jac = value_jacobian(f, x)
    return val, jac'
end

"""
    value_gradient(f::AbstractFunction, x::AbstractVector)

Returns the value and gradient of the scalar-valued function `f` at the point `x`. This is a convenience function for scalar-valued functions that simply calls [`value_jacobian_transpose`](@ref) on `f` and `x`.
"""
function value_gradient(f::AbstractFunction, x::AbstractVector)
    @assert getdim(f) == 1
    return value_jacobian_transpose(f, x)
end
