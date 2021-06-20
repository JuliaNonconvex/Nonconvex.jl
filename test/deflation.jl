import Ipopt, NLopt
using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = x[2] * sin(10x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "Deflated Ipopt" begin
    @testset "First order - $first_order" for first_order in [true, false]
        options = DeflatedOptions(
            ndeflations = 1, sub_options = IpoptOptions(),
        )
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = DeflatedAlg(IpoptAlg())
        x0 = [1.234, 2.345]
        r = Nonconvex.optimize(m, alg, x0, options = options)
        @test abs(r.solutions[1][1][2] - 2.9879) <= 1e-4
        @test abs(r.solutions[2][1][2] - 4.2435) <= 1e-4
    end
end

@testset "Deflated NLopt MMA" begin
    options = DeflatedOptions(
        ndeflations = 1, sub_options = NLoptOptions(),
    )
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    alg = DeflatedAlg(NLoptAlg(:LD_MMA))
    x0 = [1.234, 2.345]
    r = Nonconvex.optimize(m, alg, x0, options = options)
    @test abs(r.solutions[1][1][2] - 5.4996) <= 1e-4
    @test abs(r.solutions[2][1][2] - 6.7559) <= 1e-4
end
