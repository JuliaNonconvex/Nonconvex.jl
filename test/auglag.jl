using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences
const FDM = FiniteDifferences

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

options = Nonconvex.MMAOptions(tol = Nonconvex.Tolerance(kkt = 1e-6, f = 0.0))

m = Nonconvex.Model(f)
addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
add_ineq_constraint!(
    m,
    FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2),
)
x0 = [2.0, 2.0]

alg1 = MMA02()
r1 = optimize(m, alg1, x0)
@test abs(r1.minimum - sqrt(8/27)) < 1e-2
@test norm(r1.minimizer - [1/3, 8/27]) < 1e-2

alg2 = AugLag(primaloptimizer = MMA87())
r2 = optimize(m, alg2, x0)
@test abs(r2.minimum - sqrt(8/27)) < 1e-2
@test norm(r2.minimizer - [1/3, 8/27]) < 1e-2
