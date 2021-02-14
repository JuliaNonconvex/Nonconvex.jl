"""
    getasymptotes(approx::AbstractMMAApprox)

Returns the asymptotes of the MMA approximation `approx`.
"""
function getasymptotes(approx::AbstractMMAApprox)
    throw("`getasymptotes` is not defined.")
end

"""
    getapproxfg(approx::AbstractMMAApprox)

Returns the last cached value of the output of `approx`.
"""
function getapproxfg(approx::AbstractMMAApprox)
    throw("`getapproxfg` is not defined.")
end

docσ = """
    getσ(approx::AbstractMMAApprox)
    setσ!(approx::AbstractMMAApprox, σ)

Get or set the move limit `σ` of the MMA approximation `approx`. See [`MMAApprox`](@ref) for an explanation.
"""

@doc docσ
function getσ(approx::AbstractMMAApprox)
    throw("`getσ` is not defined.")
end

@doc docσ
function setσ!(approx::AbstractMMAApprox, σ)
    throw("`setσ!` is not defined.")
end

docρ = """
    getρ(approx::AbstractMMAApprox)
    setρ!(approx::AbstractMMAApprox, ρ)

Get or set the value of `ρ` in `approx`. See [`MMAApprox`](@ref) for an explanation.
"""

@doc docρ
function getρ(approx::AbstractMMAApprox)
    throw("`getρ` is not defined.")
end

@doc docρ
function setρ!(approx::AbstractMMAApprox, ρ)
    throw("`setρ!` is not defined.")
end

docfk = """
    getfk(approx::AbstractMMAApprox)
    setfk!(approx::AbstractMMAApprox, fk)

Get or set the function value `fk` of the approximated function at the approximation point `xk` of `approx`.
"""

@doc docfk
function getfk(approx::AbstractMMAApprox)
    throw("`getfk` is not defined.")
end

@doc docfk
function setfk!(approx::AbstractMMAApprox, f)
    throw("`setfk!` is not defined.")
end

docxk = """
    getxk(approx::AbstractMMAApprox)
    setxk!(approx::AbstractMMAApprox, x::AbstractVector)

Get or set the approximation point `xk` of `approx`.
"""

@doc docxk
function getxk(approx::AbstractMMAApprox)
    throw("`getfk` is not defined.")
end

@doc docxk
function setxk!(approx::AbstractMMAApprox, x::AbstractVector)
    throw("`setxk!` is not defined.")
end

doc∇fk = """
    get∇fk(approx)
    set∇fk!(approx, ∇f::AbstractVecOrMat)

Get or set the value of the gradient/Jacobian, `∇fk`, of the approximated function at the approximation point `xk`.
"""

@doc doc∇fk
function get∇fk(approx::AbstractMMAApprox)
    throw("`get∇fk` is not defined.")
end

@doc doc∇fk
function set∇fk!(approx::AbstractMMAApprox, ∇f::AbstractVecOrMat)
    throw("`set∇fk!` is not defined.")
end

docupdateapprox! = """
    updateapprox!(f, x)

Updates the MMA approximation in `f` making it around the point `x`.
"""

@doc docupdateapprox!
function updateapprox!(approx::AbstractMMAApprox, x::AbstractVector)
    throw("`updateapprox` is not defined.")
end

"""
    updateapprox!(approx::AbstractMMAApprox, x, f, ∇f)

Updates the MMA approximation `approx` making it around the point `x` where `f` is the value of the approximated function at `x` and `∇f` is its gradient if `f` is a number and the Jacobian if `f` is a vector.
"""
function updateapprox!(approx::AbstractMMAApprox, x, f, ∇f)
    throw("`updateapprox` is not defined.")
end

"""
    updateapprox!(approx::AbstractMMAApprox)

Updates the MMA approximation `approx`. This is called when `σ`, `ρ` or `xk` change.
"""
function updateapprox!(approx::AbstractMMAApprox)
    throw("`updateapprox` is not defined.")
end
