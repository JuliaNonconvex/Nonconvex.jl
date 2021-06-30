import Ipopt, NLopt, AbstractGPs
using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "Cheap objective and constraints" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 10, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Expensive objective" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        add_eq_constraint!(m, x -> sum(x) - 1/3 - 8/27)

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 10, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, skip = 2, every = 2,
            fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Expensive objective - postoptimize" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        add_eq_constraint!(m, x -> sum(x) - 1/3 - 8/27)

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 10, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = true, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-5
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-5
    end
end

@testset "Expensive objective - no preinitialize" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        add_eq_constraint!(m, x -> sum(x) - 1/3 - 8/27)

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 10, ctol = 1e-4,
            initialize = false, postoptimize = false, skip = 2, every = 2, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-5
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-5
    end
end

@testset "Expensive inequality constraint" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0), flags = [:expensive])
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 100, ctol = 1e-4,
            ninit = 5, initialize = true, postoptimize = false, last = 20, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Expensive inequality constraint with equality" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0), flags = [:expensive])
        add_ineq_constraint!(m, x -> g(x, -1, 1))

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 20, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, last = 10, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Two expensive inequality constraint" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0), flags = [:expensive])
        add_ineq_constraint!(m, x -> g(x, -1, 1), flags = [:expensive])

        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 40, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, every = 5, skip = 2,
            last = 10, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Expensive equality constraint" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        add_eq_constraint!(m, x -> sum(x) - 1/3 - 8/27, flags = [:expensive])
        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 30, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, last = 10, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(
            m, alg, [1.234, 2.345], options = options,
        )
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Expensive block constraint" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> [g(x, 2, 0), g(x, -1, 1)], flags = [:expensive])
        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 20, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, skip = 2, every = 2,
            last = 10, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(
            m, alg, [1.234, 2.345], options = options,
        )
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "All expensive no equality" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> [g(x, 2, 0), g(x, -1, 1)], flags = [:expensive])
        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(), maxiter = 50, ctol = 1e-3,
            ninit = 2, initialize = true, postoptimize = false, last = 10,
            fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "All expensive with equality" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        m = Model()
        set_objective!(m, f, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> [g(x, 2, 0), g(x, -1, 1)], flags = [:expensive])
        add_eq_constraint!(m, x -> sum(x) - 1/3 - 8/27, flags = [:expensive])
        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 50, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, last = 30, skip = 2,
            every = 2, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end

@testset "Passing surrogates in" begin
    @testset "Fit prior: $fit_prior" for fit_prior in (true, false)
        x0 = [1.234, 2.345]
        s1 = Nonconvex.surrogate(f, x0)
        s2 = Nonconvex.surrogate(x -> [g(x, 2, 0), g(x, -1, 1)], x0)

        m = Model()
        set_objective!(m, x -> s1(x).lo, flags = [:expensive])
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(
            m, x -> getproperty.(s2(x), :lo), flags = [:expensive],
        )
        alg = BayesOptAlg(IpoptAlg())
        options = BayesOptOptions(
            sub_options = IpoptOptions(print_level = 0), maxiter = 50, ctol = 1e-4,
            ninit = 2, initialize = true, postoptimize = false, fit_prior = fit_prior,
        )
        r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, surrogates = [s1, s2])
        @test abs(r.minimum - sqrt(8/27)) < 1e-3
        @test norm(r.minimizer - [1/3, 8/27]) < 1e-3
    end
end
