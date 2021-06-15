using Test, Nonconvex
import JuMP

@testset "Symbol variable names" begin
    model = JuMP.Model()
    JuMP.@variable(model, x >= 0)
    JuMP.@variable(model, 0 <= y <= 3)
    JuMP.@objective(model, Min, 20y)
    JuMP.@constraint(model, c1, 6x + 8y >= 100)
    JuMP.@constraint(model, c2, 7x + 12y >= 120)
    JuMP.@constraint(model, c3, x + y == 25)
    JuMP.@constraint(model, c4, x + y <= 50)

    dict_model = DictModel(model)
    obj = Nonconvex.getobjective(dict_model)
    ineq = Nonconvex.getineqconstraints(dict_model)
    eq = Nonconvex.geteqconstraints(dict_model)
    x = Nonconvex.getinit(dict_model)

    @test getmin(dict_model) == OrderedDict(:x => 0.0, :y => 0.0)
    @test getmax(dict_model) == OrderedDict(:x => Inf, :y => 3.0)
    @test obj(x) == 20 * x[:y]
    @test ineq(x) == [
        100 - 6 * x[:x] - 8 * x[:y],
        120 - 7 * x[:x] - 12 * x[:y],
        x[:x] + x[:y] - 50,
    ]
    @test eq(x) == [x[:x] + x[:y] - 25]

    vec_model, _, _ = Nonconvex.tovecmodel(dict_model)
    vec_x = Nonconvex.getinit(vec_model)
    vec_obj = Nonconvex.getobjective(vec_model)
    vec_ineq = Nonconvex.getineqconstraints(vec_model)
    vec_eq = Nonconvex.geteqconstraints(vec_model)

    val0, grad0 = Nonconvex.value_gradient(vec_obj, vec_x)
    val1, jac1 = Nonconvex.value_jacobian(vec_ineq, vec_x)
    val2, jac2 = Nonconvex.value_jacobian(vec_eq, vec_x)

    @test val0 == obj(x)
    @test val1 == ineq(x)
    @test val2 == eq(x)

    @test grad0 == [0, 20]
    @test jac1 == [-6 -8; -7 -12; 1 1;]
    @test jac2 == [1 1;]
end

@testset "x[ind] variable names 1" begin
    model = JuMP.Model()
    lbs = [0.0, 0.0]
    ubs = [Inf, 3.0]
    JuMP.@variable(model, x[i=1:2], lower_bound = lbs[i], upper_bound = ubs[i])
    JuMP.@objective(model, Min, 12x[1])
    JuMP.@constraint(model, c1, 6x[1] + 8x[2] >= 100)
    JuMP.@constraint(model, c2, 7x[1] + 12x[2] >= 120)
    JuMP.@constraint(model, c3, x[1] + x[2] == 25)
    JuMP.@constraint(model, c4, x[1] + x[2] <= 50)

    dict_model = DictModel(model)
    obj = Nonconvex.getobjective(dict_model)
    ineq = Nonconvex.getineqconstraints(dict_model)
    eq = Nonconvex.geteqconstraints(dict_model)
    x = Nonconvex.getinit(dict_model)

    @test getmin(dict_model) == OrderedDict(:x => [0.0, 0.0])
    @test getmax(dict_model) == OrderedDict(:x => [Inf, 3.0])
    @test obj(x) == 12 * x[:x][1]
    @test ineq(x) == [
        100 - 6 * x[:x][1] - 8 * x[:x][2],
        120 - 7 * x[:x][1] - 12 * x[:x][2],
        x[:x][1] + x[:x][2] - 50,
    ]
    @test eq(x) == [x[:x][1] + x[:x][2] - 25]

    vec_model, _, _ = Nonconvex.tovecmodel(dict_model)
    vec_x = Nonconvex.getinit(vec_model)
    vec_obj = Nonconvex.getobjective(vec_model)
    vec_ineq = Nonconvex.getineqconstraints(vec_model)
    vec_eq = Nonconvex.geteqconstraints(vec_model)

    val0, grad0 = Nonconvex.value_gradient(vec_obj, vec_x)
    val1, jac1 = Nonconvex.value_jacobian(vec_ineq, vec_x)
    val2, jac2 = Nonconvex.value_jacobian(vec_eq, vec_x)

    @test val0 == obj(x)
    @test val1 == ineq(x)
    @test val2 == eq(x)

    @test grad0 == [12, 0]
    @test jac1 == [-6 -8; -7 -12; 1 1;]
    @test jac2 == [1 1;]
end

@testset "x[ind] variable names 2" begin
    model = JuMP.Model()
    lbs = [0.0, 0.0]
    ubs = [Inf, 3.0]
    JuMP.@variable(model, x[i=2:-1:1], lower_bound = lbs[i], upper_bound = ubs[i])
    JuMP.@objective(model, Min, 12x[1])
    JuMP.@constraint(model, c1, 6x[1] + 8x[2] >= 100)
    JuMP.@constraint(model, c2, 7x[1] + 12x[2] >= 120)
    JuMP.@constraint(model, c3, x[1] + x[2] == 25)
    JuMP.@constraint(model, c4, x[1] + x[2] <= 50)

    dict_model = DictModel(model)
    obj = Nonconvex.getobjective(dict_model)
    ineq = Nonconvex.getineqconstraints(dict_model)
    eq = Nonconvex.geteqconstraints(dict_model)
    x = Nonconvex.getinit(dict_model)

    @test getmin(dict_model) == OrderedDict(:x => JuMP.Containers.DenseAxisArray(lbs[end:-1:1], 2:-1:1))
    @test getmax(dict_model) == OrderedDict(:x => JuMP.Containers.DenseAxisArray(ubs[end:-1:1], 2:-1:1))
    @test obj(x) == 12 * x[:x][1]
    @test ineq(x) == [
        100 - 6 * x[:x][1] - 8 * x[:x][2],
        120 - 7 * x[:x][1] - 12 * x[:x][2],
        x[:x][1] + x[:x][2] - 50,
    ]
    @test eq(x) == [x[:x][1] + x[:x][2] - 25]

    vec_model, _, _ = Nonconvex.tovecmodel(dict_model)
    vec_x = Nonconvex.getinit(vec_model)
    vec_obj = Nonconvex.getobjective(vec_model)
    vec_ineq = Nonconvex.getineqconstraints(vec_model)
    vec_eq = Nonconvex.geteqconstraints(vec_model)

    val0, grad0 = Nonconvex.value_gradient(vec_obj, vec_x)
    val1, jac1 = Nonconvex.value_jacobian(vec_ineq, vec_x)
    val2, jac2 = Nonconvex.value_jacobian(vec_eq, vec_x)

    @test val0 == obj(x)
    @test val1 == ineq(x)
    @test val2 == eq(x)

    @test grad0 == [0, 12]
    @test jac1 == [-8 -6; -12 -7; 1 1;]
    @test jac2 == [1 1;]
end

@testset "x[(ind1, ind2)] variable names" begin
    model = JuMP.Model()
    lbs = [0.0, 0.0]
    ubs = [Inf, 3.0]
    inds = [(1,1), (2,1)]
    JuMP.@variable(model, x[ind = inds], lower_bound = lbs[ind[1]], upper_bound = ubs[ind[1]])
    JuMP.@objective(model, Min, 12x[(1,1)] + 20x[(2,1)])
    JuMP.@constraint(model, c1, 6x[(1,1)] + 8x[(2,1)] >= 100)
    JuMP.@constraint(model, c2, 7x[(1,1)] + 12x[(2,1)] >= 120)
    JuMP.@constraint(model, c3, x[(1,1)] + x[(2,1)] == 25)
    JuMP.@constraint(model, c4, x[(1,1)] + x[(2,1)] <= 50)

    dict_model = DictModel(model)
    obj = Nonconvex.getobjective(dict_model)
    ineq = Nonconvex.getineqconstraints(dict_model)
    eq = Nonconvex.geteqconstraints(dict_model)
    x = Nonconvex.getinit(dict_model)

    @test getmin(dict_model) == OrderedDict(:x => JuMP.Containers.DenseAxisArray(lbs, [(1, 1), (2, 1)]))
    @test getmax(dict_model) == OrderedDict(:x => JuMP.Containers.DenseAxisArray(ubs, [(1, 1), (2, 1)]))
    @test obj(x) == 12 * x[:x][(1, 1)] + 20 * x[:x][(2, 1)]
    @test ineq(x) == [
        100 - 6 * x[:x][(1, 1)] - 8 * x[:x][(2, 1)],
        120 - 7 * x[:x][(1, 1)] - 12 * x[:x][(2, 1)],
        x[:x][(1, 1)] + x[:x][(2, 1)] - 50,
    ]
    @test eq(x) == [x[:x][(1, 1)] + x[:x][(2, 1)] - 25]

    vec_model, _, _ = Nonconvex.tovecmodel(dict_model)
    vec_x = Nonconvex.getinit(vec_model)
    vec_obj = Nonconvex.getobjective(vec_model)
    vec_ineq = Nonconvex.getineqconstraints(vec_model)
    vec_eq = Nonconvex.geteqconstraints(vec_model)

    val0, grad0 = Nonconvex.value_gradient(vec_obj, vec_x)
    val1, jac1 = Nonconvex.value_jacobian(vec_ineq, vec_x)
    val2, jac2 = Nonconvex.value_jacobian(vec_eq, vec_x)

    @test val0 == obj(x)
    @test val1 == ineq(x)
    @test val2 == eq(x)

    @test grad0 == [12, 20]
    @test jac1 == [-6 -8; -7 -12; 1 1;]
    @test jac2 == [1 1;]
end
