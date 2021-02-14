using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences
const FDM = FiniteDifferences

f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a*x[1] + b)^3 - x[2]

m = Model(f)
addvar!(m, [0.0, 0.0], [10.0, 10.0])
add_ineq_constraint!(m, x -> g(x, 2, 0))
add_ineq_constraint!(m, x -> g(x, -1, 1))

@test getmin(m) == fill(0, 2)
@test getmax(m) == fill(10, 2)
@test getmin(m, 1) == 0.0
@test getmax(m, 1) == 10.0
@test getmin(m, 2) == 0.0
@test getmax(m, 2) == 10.0

p = [1.234, 2.345]

val0, grad0 = Nonconvex.value_gradient(Nonconvex.getobjective(m), p)
@test val0 == f(p)
@test grad0 == Zygote.gradient(f, p)[1]

val1, grad1 = Nonconvex.value_gradient(Nonconvex.getineqconstraint(m, 1), p)
@test val1 == g(p, 2, 0)
@test grad1 == Zygote.gradient(p -> g(p, 2, 0), p)[1]

val2, grad2 = Nonconvex.value_gradient(Nonconvex.getineqconstraint(m, 2), p)
@test val2 == g(p, -1, 1)
@test grad2 == Zygote.gradient(p -> g(p, -1, 1), p)[1]

vals, jac = Nonconvex.value_jacobian(Nonconvex.getineqconstraints(m), p)
@test [val1, val2] == vals
@test [grad1 grad2]' == jac

vals, jac = Nonconvex.value_jacobian(Nonconvex.getobjectiveconstraints(m), p)
@test [val0, val1, val2] == vals
@test [grad0 grad1 grad2]' == jac
