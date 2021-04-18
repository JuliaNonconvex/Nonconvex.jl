using Nonconvex, Hyperopt
using SafeTestsets, Test

f(x) = sqrt(x[2])
g(x, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0])
add_ineq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))

alg = MMA87() # or MMA02()
options = Nonconvex.MMAOptions(
    tol = Nonconvex.Tolerance(kkt = 1e-6, f = 0.0), s_init = 0.1,
)
convcriteria = KKTCriteria()

function is_best_result(target_result::Float64, results::Array{Nonconvex.AbstractResult})
    for _result in results
        if _result.minimum < target_result
            false
        end
    end
    true
end

@testset "search_x0" begin
    r1 = @search_x0 Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
    r2 = @search_x0 X0OptOptions(), Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
    r3 = @hypersearch X0OptOptions(), Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
    
    # Customized options
    hyperopt_options = X0OptOptions(x0_lb=[0.5, 0.5], x0_rb=[2.8, 2.8],
    searchspace_size=1000, iters=20, 
    sampler=Hyperopt.RandomSampler(), 
    verbose=true,
    keepall=true)
    # Searching hyperparameters using customized options. 
    r4 = @hypersearch hyperopt_options, Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
    # Equivalent as above. 
    r5 = @search_x0 hyperopt_options, Nonconvex.optimize(m, alg, [1.234, 2.345], options = options, convcriteria = convcriteria)
    for r in [r1, r2, r3, r4, r5]
        @test is_best_result(r.minimum, r.results)
    end
end
