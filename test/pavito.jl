import Pavito
using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "First order - $first_order" for first_order in [true, false]
    options = PavitoIpoptCbcOptions(first_order = first_order)
    @testset "Simple constraints" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = PavitoIpoptCbcAlg()

        setinteger!(m, 1, true)
        r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r1.minimizer[1] - round(Int, r1.minimizer[1]) ≈ 0 atol = 1e-7

        setinteger!(m, 1, false)
        setinteger!(m, 2, true)
        r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r2.minimizer[2] - round(Int, r2.minimizer[2]) ≈ 0 atol = 1e-7
    end

    @testset "Block constraints" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2))

        alg = PavitoIpoptCbcAlg()

        setinteger!(m, 1, true)
        r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r1.minimizer[1] - round(Int, r1.minimizer[1]) ≈ 0 atol = 1e-7

        setinteger!(m, 1, false)
        setinteger!(m, 2, true)
        r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r2.minimizer[2] - round(Int, r2.minimizer[2]) ≈ 0 atol = 1e-7
    end

    @testset "Infinite bounds" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = PavitoIpoptCbcAlg()

        setinteger!(m, 1, true)
        r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r1.minimizer[1] - round(Int, r1.minimizer[1]) ≈ 0 atol = 1e-7

        setinteger!(m, 1, false)
        setinteger!(m, 2, true)
        r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r2.minimizer[2] - round(Int, r2.minimizer[2]) ≈ 0 atol = 1e-7
    end
end
