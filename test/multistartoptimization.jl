using Nonconvex, LinearAlgebra, Test

f(x) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

@testset "Simple constraints" begin
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    mso_optimizer = MsOAlg(TikTak(100))
    mso_local_method = NLopt.LN_BOBYQA
    mso_options = MsOOptions(use_threads = false)
    mso_res = get_starting_point(m, mso_optimizer, local_method = mso_local_method, options = mso_options)

    @test mso_res.minimum == 0.0
end

@testset "Without upper or lower bounds" begin
    try
        m = Model(f)
        mso_optimizer = MsOAlg(TikTak(100))
        mso_local_method = NLopt.LN_BOBYQA
        mso_res = get_starting_point(m, mso_optimizer, local_method = mso_local_method)
    catch err
        @test isa(err, ErrorException)
    end
end
