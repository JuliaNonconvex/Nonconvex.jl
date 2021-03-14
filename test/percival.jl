using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

options = PercivalOptions()

@testset "Simple constraints" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    alg = PercivalAlg()
    r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
    @test abs(r.minimum - sqrt(8/27)) < 1e-6
    @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
end

@testset "Block constraints" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2))

    alg = PercivalAlg()
    r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
    @test abs(r.minimum - sqrt(8/27)) < 1e-6
    @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
end

@testset "Infinite bounds" begin
    @testset "Infinite upper bound" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = PercivalAlg()
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-6
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
    end
    @testset "Infinite lower bound" begin
        m = Model(f)
        addvar!(m, [-Inf, -Inf], [10, 10])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = PercivalAlg()
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-6
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
    end
    @testset "Infinite upper and lower bound" begin
        m = Model(f)
        addvar!(m, [-Inf, -Inf], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = PercivalAlg()
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-6
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
    end
end
