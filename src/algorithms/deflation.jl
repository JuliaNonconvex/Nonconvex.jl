# The DeflationOperator implementaion is adapted from BifurcationKit.jl which has the following license:
#=
The BifurcationKit.jl Julia package is licensed under the MIT "Expat" License:
    Copyright (c) 2019: Romain VELTZ.
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
struct DeflationOperator{T <: Real, Tdot, vectype}
    radius::T
	power::T
	dot::Tdot
	shift::T
	roots::Vector{vectype}
end
function DeflationOperator(
    X; radius = 0.0, power = 2.0,
    dot = LinearAlgebra.dot, shift = 1.0,
)
    return DeflationOperator(radius, power, dot, shift, X)
end
function (df::DeflationOperator{T, Tdot, vectype})(u::vectype) where {T, Tdot, vectype}
	nrm  = u -> df.dot(u, u)
	if length(df.roots) == 0
		return T(1)
	end
	# compute u - df.roots[1]
	tmp = 1 * u
    axpy!(T(-1), df.roots[1], tmp)
	out = T(1) / max(nrm(tmp) - df.radius, 0)^df.power + df.shift
	for ii in 2:length(df.roots)
		copyto!(tmp, u)
        axpy!(T(-1), df.roots[ii], tmp)
		out *= T(1) / max(nrm(tmp) - df.radius, 0)^df.power + df.shift
	end
	return out
end
function ChainRulesCore.rrule(f::DeflationOperator, x)
    return f.def(x), Δ -> (nothing, f.def(x) * pb(Δ)[1])
end

@params struct DeflatedFunction <: Function
    f
    def
end
(f::DeflatedFunction)(x) = f.f(x)
function ChainRulesCore.rrule(f::DeflatedFunction, x)
    v, pb = Zygote.pullback(f.f, x)
    return v, Δ -> (nothing, f.def(x) * pb(Δ)[1])
end

@params struct DeflatedOptions
    ndeflations::Int
    radius::Float64
	power::Float64
	dot::Function
	shift::Float64
    sub_options
end
function DeflatedOptions(;
    ndeflations = 2,
    radius = 0.0,
	power = 2.0,
	dot = LinearAlgebra.dot,
	shift = 1.0,
    sub_options = IpoptOptions(),
    kwargs...,
)
    return DeflatedOptions(
        ndeflations, radius, power, dot, shift, sub_options,
    )
end

struct DeflatedAlg{A} <: AbstractOptimizer
    sub_alg::A
end

function deflatedmodel(vecmodel::VecModel, deflate)
    _obj = vecmodel.objective
    if _obj isa Function
        obj = Objective(
            DeflatedFunction(_obj.f, deflate),
            _obj.multiple,
            _obj.flags,
        )
    else
        obj = _obj
    end
    ineq_constraints = map(vecmodel.ineq_constraints.fs) do c
        @assert c isa IneqConstraint
        return IneqConstraint(
            DeflatedFunction(c, deflate), zero.(c.rhs), c.dim, c.flags,
        )
    end |> VectorOfFunctions
    eq_constraints = map(vecmodel.eq_constraints.fs) do c
        @assert c isa EqConstraint
        return EqConstraint(
            DeflatedFunction(c, deflate), zero.(c.rhs), c.dim, c.flags,
        )
    end |> VectorOfFunctions
    return VecModel(
        obj, eq_constraints, ineq_constraints, vecmodel.box_min,
        vecmodel.box_max, vecmodel.init, vecmodel.integer,
    )
end

@params mutable struct DeflatedWorkspace <: Workspace
    model::VecModel
    defmodel::VecModel
    defop::DeflationOperator
    sub_workspace
    x0::AbstractVector
    X::AbstractVector
    options::DeflatedOptions
end
function Workspace(
    model::VecModel, optimizer::DeflatedAlg, x0::AbstractVector;
    options::DeflatedOptions = DeflatedOptions(), kwargs...,
)
    @unpack radius, shift, power, dot, sub_options = options
    X = Vector{Float64}[]
    defop = DeflationOperator(
        X, radius = radius, shift = shift,
        dot = dot, power = power,
    )
    defmodel = deflatedmodel(model, defop)
    sub_workspace = Workspace(
        defmodel, optimizer.sub_alg, copy(x0);
        options = sub_options, kwargs...,
    )
    return DeflatedWorkspace(
        model, defmodel, defop, sub_workspace, x0, X, options,
    )
end

@params struct DeflatedResult
    solutions
    distances
    minimizer
    minimum
end

function optimize!(workspace::DeflatedWorkspace)
    @unpack options, x0, sub_workspace = workspace
    @unpack model, defop, X = workspace
    solutions = Tuple{Vector{Float64}, Float64}[]
    @unpack ndeflations = options
    sub_workspace.x0 .= x0
    r = optimize!(sub_workspace)
    minimizer = copy(r.minimizer)
    minval = getobjective(model)(minimizer)
    push!(X, copy(minimizer))
    push!(solutions, (copy(minimizer), minval))
    if ndeflations > 0
        i = 0
        while i < ndeflations
            sub_workspace.x0 .= x0
            r = optimize!(sub_workspace)
            objval = getobjective(model)(r.minimizer)
            push!(X, copy(r.minimizer))
            push!(solutions, (copy(r.minimizer), objval))
            if objval <= minval
                minval = objval
                minimizer = copy(r.minimizer)
            end
            i += 1
        end
    end
    distances = [
        norm(solutions[i][1] - solutions[j][1]) for
        i in 1:(ndeflations + 1), j in 1:(ndeflations + 1)
    ]
    return DeflatedResult(
        solutions, distances, minimizer, minval,
    )
end
