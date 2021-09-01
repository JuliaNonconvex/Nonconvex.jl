using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences
const FDM = FiniteDifferences

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

options = MMAOptions(
    tol = Tolerance(kkt = 1e-6, f = 0.0),
    s_init = 0.1,
)

@testset "Simple constraints" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    @testset "MMA $(alg isa MMA87 ? "1987" : "2002")" for alg in (MMA87(), MMA02())
        for convcriteria in (KKTCriteria(), IpoptCriteria())
            r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
            @test abs(r.minimum - sqrt(8/27)) < 1e-6
            @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
        end
    end
end

@testset "Block constraints" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2))

    @testset "MMA $(alg isa MMA87 ? "1987" : "2002")" for alg in (MMA87(), MMA02())
        for convcriteria in (KKTCriteria(), IpoptCriteria())
            r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
            @test abs(r.minimum - sqrt(8/27)) < 1e-6
            @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
        end
    end
end

@testset "Infinite bounds" begin
    @testset "Infinite upper bound" begin
        m = Model(f)
        addvar!(m, [0.0, 0.0], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        @testset "MMA $(alg isa MMA87 ? "1987" : "2002")" for alg in (MMA87(), MMA02())
            for convcriteria in (KKTCriteria(), IpoptCriteria())
                r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
                @test abs(r.minimum - sqrt(8/27)) < 1e-6
                @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
            end
        end
    end
    @testset "Infinite lower bound" begin
        m = Model(f)
        addvar!(m, [-Inf, -Inf], [10, 10])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        @testset "MMA $(alg isa MMA87 ? "1987" : "2002")" for alg in (MMA87(), MMA02())
            for convcriteria in (KKTCriteria(), IpoptCriteria())
                r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
                @test abs(r.minimum - sqrt(8/27)) < 1e-6
                @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
            end
        end
    end
    @testset "Infinite upper and lower bound" begin
        m = Model(f)
        addvar!(m, [-Inf, -Inf], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        @testset "MMA $(alg isa MMA87 ? "1987" : "2002")" for alg in (MMA87(), MMA02())
            for convcriteria in (KKTCriteria(), IpoptCriteria())
                r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
                @test abs(r.minimum - sqrt(8/27)) < 1e-6
                @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
            end
        end
    end
end

# Extended formulation
# Not great for non-convex problems so we need to start pretty close
# The problem is that an optimum solution to the extended formulation may not be feasible in the original formulation. This could be a local optimum or just due to a subspace of possible optimal solutions along a subspace.
@testset "Extended formulation" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    for convcriteria in (KKTCriteria(), IpoptCriteria())
        r = Nonconvex.optimize(m, MMA02(), [0.4, 0.5], options = options, convcriteria = convcriteria)
        @test abs(r.minimum - sqrt(8/27)) < 1e-6
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
    end
end

# Test unflatten minimizer and initial_x
@testset "Unflatten initial_x" begin
    function _f((x, _))
        sqrt(x[2])
    end
    function _g((x, _), a, b)
        (a*x[1] + b)^3 - x[2]
    end

    _options = MMAOptions(
        tol = Tolerance(kkt = 1e-6, f = 0.0),
        s_init = 0.1,
        maxiter = 1
    )

    m = Model(_f)
    addvar!(m, [(0.0, 0.0), (0.0, 0.0)], [(10.0, Inf), (10.0, Inf)])

    function eqc1(x)
        _g(x, 2, 0)
    end
    function eqc2(x)
        _g(x, -1, 1)
    end
    add_ineq_constraint!(m, eqc1)
    add_ineq_constraint!(m, eqc2)

    x0 = [(0.0, 0.0), (0.0, 0.0)]
    r = Nonconvex.optimize(m, MMA02(), x0, options = _options, convcriteria = KKTCriteria())
    @test size(r.minimizer) == (2,)
    @test size(r.initial_x) == (2,)
end
