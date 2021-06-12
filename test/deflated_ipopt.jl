import Ipopt
using Nonconvex, LinearAlgebra, Test

#f(x::AbstractVector) = x[2] * sin(10x[2])
_f(x) = (x - 0.1) * (x - 0.3) * (x - 0.5) * (x - 1.0)
f(x::AbstractVector) = _f(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "First order - $first_order" for first_order in [true, false]
    options = DeflatedIpoptOptions(first_order = first_order, niters = 5, tol = 1e-4)
    @testset "Simple constraints" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = DeflatedIpoptAlg()
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
    end
end
