@params struct DeflatedIpoptOptions
    niters::Int
    ipopt_options::IpoptOptions
end
function DeflatedIpoptOptions(;
    niters = 1,
    first_order = true,
    hessian_approximation = first_order ? "limited-memory" : "exact",
    kwargs...,
)
    h = hessian_approximation
    return DeflatedIpoptOptions(
        niters,
        IpoptOptions((;hessian_approximation = h, kwargs...)),
    )
end

@params mutable struct DeflatedIpoptWorkspace <: Workspace
    model::VecModel
    problem::Ipopt.IpoptProblem
    x0::AbstractVector
    X::AbstractVector
    options::DeflatedIpoptOptions
    counter::Base.RefValue{Int}
end
function DeflatedIpoptWorkspace(
    model::VecModel, x0::AbstractVector = getinit(model), X = Vector{Float64}[];
    options = DeflatedIpoptOptions(), kwargs...,
)
    problem, counter = get_deflated_ipopt_problem(
        model, copy(x0),
        options.ipopt_options.nt.hessian_approximation == "limited-memory",
        X,
    )
    return DeflatedIpoptWorkspace(model, problem, x0, X, options, counter)
end
@params struct DeflatedIpoptResult
    solutions
    distances
    minimizer
    minimum
    problem
    status
    fcalls::Int
end

function optimize!(workspace::DeflatedIpoptWorkspace)
    @unpack problem, options, counter, x0, X, model = workspace
    solutions = Tuple{Vector{Float64}, Float64}[]
    niters = options.niters
    counter[] = 0
    foreach(keys(options.ipopt_options.nt)) do k
        v = options.ipopt_options.nt[k]
        Ipopt.addOption(problem, string(k), v)
    end
    solvestat = Ipopt.solveProblem(problem)
    minimizer = copy(problem.x)
    minval = getobjective(model)(problem.x)
    status = solvestat
    push!(X, copy(problem.x))
    push!(solutions, (copy(problem.x), minval))
    if niters > 1
        i = 1
        while i < niters
            problem.x .= x0
            solvestat = Ipopt.solveProblem(problem)
            objval = getobjective(model)(problem.x)
            #if solvestat == -1
            #    break
            #=else=#if objval < minval
                status = solvestat
                minval = objval
                minimizer = copy(problem.x)
            end
            push!(X, copy(problem.x))
            push!(solutions, (copy(problem.x), objval))
            i += 1
        end
    end
    distances = [
        norm(solutions[i][1] - solutions[j][1]) for
        i in 1:niters, j in 1:niters
    ]
    return DeflatedIpoptResult(
        solutions, distances, minimizer, minval, problem, status, counter[],
    )
end

struct DeflatedIpoptAlg{O} <: AbstractOptimizer
    options::O
end
DeflatedIpoptAlg(; kwargs...) = DeflatedIpoptAlg(kwargs)

function Workspace(model::VecModel, optimizer::DeflatedIpoptAlg, args...; kwargs...,)
    return DeflatedIpoptWorkspace(model, args...; kwargs...)
end

# The DeflationOperator implementaion is adapted from BifurcationKit.jl which has the following license:
#=
The BifurcationKit.jl Julia package is licensed under the MIT "Expat" License:
    Copyright (c) 2019: Romain VELTZ.
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
struct DeflationOperator{T <: Real, Tdot, vectype}
	power::T
	dot::Tdot
	shift::T
	roots::Vector{vectype}
end
function DeflationOperator(X)
    return DeflationOperator(2.0, dot, 1.0, X)
end

function (df::DeflationOperator{T, Tdot, vectype})(u::vectype) where {T, Tdot, vectype}
    r = 10.0
	nrm  = u -> df.dot(u, u)
	if length(df.roots) == 0
		return T(1)
	end
	# compute u - df.roots[1]
	tmp = _copy(u);	axpy!(T(-1), df.roots[1], tmp)
	out = T(1) / max(nrm(tmp) - r, 0)^df.power + df.shift
	for ii in 2:length(df.roots)
		copyto!(tmp, u); axpy!(T(-1), df.roots[ii], tmp)
		out *= T(1) / max(nrm(tmp) - r, 0)^df.power + df.shift
	end
	return out
end
_copy(b) = 1*b

function get_deflated_ipopt_problem(model::VecModel, x0::AbstractVector, first_order::Bool, X)
    eq = if length(model.eq_constraints.fs) == 0
        nothing
    else
        model.eq_constraints
    end
    ineq = if length(model.ineq_constraints.fs) == 0
        nothing
    else
        model.ineq_constraints
    end
    obj = CountingFunction(getobjective(model))
    return get_deflated_ipopt_problem(
        obj,
        ineq,
        eq,
        x0,
        getmin(model),
        getmax(model),
        first_order,
        X,
    ), obj.counter
end
function get_deflated_ipopt_problem(obj, ineq_constr, eq_constr, x0, xlb, xub, first_order, X)
    deflation_op = DeflationOperator(X)
    nvars = 0
    if ineq_constr !== nothing
        ineqJ0 = Zygote.jacobian(ineq_constr, x0)[1]
        ineq_nconstr, nvars = size(ineqJ0)
        Joffset = nvalues(ineqJ0)
    else
        ineqJ0 = nothing
        ineq_nconstr = 0
        Joffset = 0
    end
    if eq_constr !== nothing
        eqJ0 = Zygote.jacobian(eq_constr, x0)[1]
        eq_nconstr, nvars = size(eqJ0)
    else
        eqJ0 = nothing
        eq_nconstr = 0
    end
    @assert nvars > 0
    lag(factor, y) = x -> begin
        factor * obj(x) + 
            _dot(ineq_constr, x, @view(y[1:ineq_nconstr])) + 
            _dot(eq_constr, x, @view(y[ineq_nconstr+1:end]))
    end
    clb = [fill(-Inf, ineq_nconstr); zeros(eq_nconstr)]
    cub = zeros(ineq_nconstr + eq_nconstr)

    function eval_g(x::Vector{Float64}, g::Vector{Float64})
        if ineq_constr !== nothing
            g[1:ineq_nconstr] .= ineq_constr(x)
        end
        if eq_constr !== nothing
            g[ineq_nconstr+1:end] .= eq_constr(x)
        end
        return g
    end
    function eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
        grad_f .= deflation_op(x) .* Zygote.gradient(obj, x)[1]
    end
    function eval_jac_g(x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector{Float64})
        if mode == :Structure
            ineqJ0 === nothing || fill_indices!(rows, cols, ineqJ0)
            eqJ0 === nothing || fill_indices!(rows, cols, eqJ0, offset = Joffset)
        else
            values .= 0
            if ineq_constr !== nothing
                ineqJ = deflation_op(x) .* Zygote.jacobian(ineq_constr, x)[1]
                add_values!(values, ineqJ)
            end
            if eq_constr !== nothing
                eqJ = deflation_op(x) .* Zygote.jacobian(eq_constr, x)[1]
                add_values!(values, eqJ, offset = Joffset)
            end
        end
    end

    if first_order
        eval_h = (x...) -> 0.0
        Hnvalues = 0
    else
        HL0 = LowerTriangular(
            Zygote.hessian(
                lag(1.0, ones(ineq_nconstr + eq_nconstr)),
                x0,
            ),
        )
        eval_h = function (x::Vector{Float64}, mode, rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Vector{Float64})
            if mode == :Structure
                fill_indices!(rows, cols, HL0)
            else
                HL = LowerTriangular(
                    deflation_op(x)^2 .* Zygote.hessian(lag(obj_factor, lambda), x),
                )
                values .= 0
                add_values!(values, HL)
            end
        end
        Hnvalues = nvalues(HL0)
    end
    prob = Ipopt.createProblem(
        nvars, xlb, xub, ineq_nconstr + eq_nconstr, clb, cub,
        nvalues(ineqJ0) + nvalues(eqJ0), Hnvalues, obj,
        eval_g, eval_grad_f, eval_jac_g, eval_h,
    )
    prob.x = x0
    return prob
end
