import Juniper
using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "First order - $first_order" for first_order in [true, false]
    options = JuniperIpoptOptions(first_order = first_order)
    @testset "Simple constraints" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = JuniperIpoptAlg()
        r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r1.minimum - sqrt(8/27)) < 1e-6
        @test norm(r1.minimizer - [1/3, 8/27]) < 1e-6

        setinteger!(m, 1, true)
        r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r2.minimizer[1] - round(Int, r2.minimizer[1]) ≈ 0 atol = 1e-7

        setinteger!(m, 1, false)
        setinteger!(m, 2, true)
        r3 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r3.minimizer[2] - round(Int, r3.minimizer[2]) ≈ 0 atol = 1e-7
    end

    @testset "Block constraints" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [10.0, 10.0])
        add_ineq_constraint!(m, FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2))

        alg = JuniperIpoptAlg()
        r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r1.minimum - sqrt(8/27)) < 1e-6
        @test norm(r1.minimizer - [1/3, 8/27]) < 1e-6

        setinteger!(m, 1, true)
        r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r2.minimizer[1] - round(Int, r2.minimizer[1]) ≈ 0 atol = 1e-7

        setinteger!(m, 1, false)
        setinteger!(m, 2, true)
        r3 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test r3.minimizer[2] - round(Int, r3.minimizer[2]) ≈ 0 atol = 1e-7
    end

    @testset "Infinite bounds" begin
        @testset "Infinite upper bound" begin
            m = Model(f)
            addvar!(m, [0.0, 0.0], [Inf, Inf])
            add_ineq_constraint!(m, x -> g(x, 2, 0))
            add_ineq_constraint!(m, x -> g(x, -1, 1))

            alg = JuniperIpoptAlg()
            r1 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
            @test abs(r1.minimum - sqrt(8/27)) < 1e-6
            @test norm(r1.minimizer - [1/3, 8/27]) < 1e-6
    
            setinteger!(m, 1, true)
            r2 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
            @test r2.minimizer[1] - round(Int, r2.minimizer[1]) ≈ 0 atol = 1e-7
    
            setinteger!(m, 1, false)
            setinteger!(m, 2, true)
            r3 = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
            @test r3.minimizer[2] - round(Int, r3.minimizer[2]) ≈ 0 atol = 1e-7
        end
        @testset "Infinite lower bound" begin
            m = Model(f)
            addvar!(m, [-Inf, -Inf], [10, 10])
            add_ineq_constraint!(m, x -> g(x, 2, 0))
            add_ineq_constraint!(m, x -> g(x, -1, 1))

            alg = JuniperIpoptAlg()
            r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
            @test abs(r.minimum - sqrt(8/27)) < 1e-6
            @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
        end
        @testset "Infinite upper and lower bound" begin
            m = Model(f)
            addvar!(m, [-Inf, -Inf], [Inf, Inf])
            add_ineq_constraint!(m, x -> g(x, 2, 0))
            add_ineq_constraint!(m, x -> g(x, -1, 1))

            alg = JuniperIpoptAlg()
            r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
            @test abs(r.minimum - sqrt(8/27)) < 1e-6
            @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
        end
    end
end
