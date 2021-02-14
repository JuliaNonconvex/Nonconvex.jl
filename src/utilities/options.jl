"""
    MMAOptions

A struct that stores all the options of the MMA algorithms. Th following are the fields of `MMAOptions`:
 - `maxiter`: the maximum number of inner iterations. For `MMA87`, there is 1 inner iteration per outer iteration.
 - `outer_maxiter`: the maximum number of outer iterations.
 - `maxinner`: the maximum number of inner iterations per outer iteration of [`MMA02`](@ref). Not applicable for [`MMA87`](@ref).
 - `tol`: a tolerance struct of type [`Tolerance`](@ref).
 - `s_init`: defined in the original [`MMA02`](@ref) paper.
 - `s_incr`: defined in the original [`MMA02`](@ref) paper.
 - `s_decr`: defined in the original [`MMA02`](@ref) paper.
 - `store_trace`: if true, a trace will be stored.
 - `dual_options`: the options passed to the dual optimizer from [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl).
"""
@with_kw mutable struct MMAOptions{T, Ttol <: Tolerance, TSubOptions <: Optim.Options}
    maxiter::Int = 1000
    outer_maxiter::Int = 10^8
    maxinner::Int = 10
    tol::Ttol = Tolerance()
    s_init::T = 0.5
    s_incr::T = 1.2
    s_decr::T = 0.7
    store_trace::Bool = false
    show_trace::Bool = false
    auto_scale::Bool = false
    keep_best::Bool = false
    dual_options::TSubOptions = Optim.Options(allow_f_increases = false, iterations = 1000, outer_iterations=1000)
end
