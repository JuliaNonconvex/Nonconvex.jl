import Ipopt, Hyperopt
using Nonconvex, LinearAlgebra, Test

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

alg = HyperoptAlg(IpoptAlg())

@testset "Sampler - $spl_name" for (spl_name, spl) in [
    # ("RandomSampler", RandomSampler()),
    ("Hyperband", Hyperband(R=100, η=3, inner=RandomSampler())),
    ("LHSampler", LHSampler()),
    ("CLHSampler", CLHSampler()),
    # ("GPSampler", GPSampler()),
]
    if spl_name == "Hyperband"
        options = HyperoptOptions(
            sub_options = max_iter -> IpoptOptions(first_order = true, max_iter = max_iter),
            sampler = spl,
        )
    else
        options = HyperoptOptions(
            sub_options = IpoptOptions(first_order = true),
            sampler = spl,
        )
    end
    m = Model(f)
    addvar!(m, [0.0, 0.0], [10.0, 10.0])
    add_ineq_constraint!(m, x -> g(x, 2, 0))
    add_ineq_constraint!(m, x -> g(x, -1, 1))

    r = Nonconvex.optimize(m, alg, [1.234, 2.345], options = options)
    @test abs(r.minimum - sqrt(8/27)) < 1e-6
    @test norm(r.minimizer - [1/3, 8/27]) < 1e-6
end
