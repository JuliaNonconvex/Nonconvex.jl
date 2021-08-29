# Test cases of semidefinite programming
using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences, Distributions
using ChainRulesTestUtils, Random
Random.seed!(1)

const FDM = FiniteDifferences

# Test setting
mat_dim = 2
mat_length = mat_dim^2
n_sample = 300

@testset "Test matrix to vectors function" begin
    test_rrule(Nonconvex.rearrange_x, rand((mat_dim^2-mat_dim)÷2), rand(mat_dim))
end

function random_psd_mat(mat_dim)
    _mat = randn(mat_dim, mat_dim)
    return _mat' * _mat
end
function sd_constraint((x_L, x_D))
    decompress_symmetric(x_L, x_D)
end

@testset "Single semi-definite constraint" begin
    # Randomize groundtruth
    μ, Σ = randn(mat_dim), random_psd_mat(mat_dim)

    # Generate samples
    samples = rand(MvNormal(μ, Σ), n_sample)

    # Objective function
    function f((x_L, x_D))
        -loglikelihood(MvNormal(μ, decompress_symmetric(x_L, x_D)), samples)
    end

    model = Model(f)

    mat_x0 = random_psd_mat(mat_dim)
    x0 = [mat_x0[Nonconvex.lowertriangind(mat_x0)], diag(mat_x0)]
    lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
    ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
    addvar!(model, lbs, ubs)

    add_sd_constraint!(model, sd_constraint)

    alg = SDPBarrierAlg(sub_alg=IpoptAlg())
    options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
    result = optimize(model, alg, x0, options = options)

    minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
    _Σ = sd_constraint(minimizer)

    println("result: \n $result")

    println("minimum: \n $minimum")
    println("minimizer: \n $minimizer")
    println("_Σ: \n $(_Σ)")

    println("Σ: \n $Σ")
    println("abs(_Σ - Σ): \n $(abs.(_Σ - Σ))")
    println("mean(abs(_Σ - Σ)): \n $(mean(abs.(_Σ - Σ)))")

    @test Σ ≈ _Σ rtol = 0.1
end

@testset "Two semi-definite constraints" begin
    # Randomize groundtruth
    μ1, Σ1 = randn(mat_dim), random_psd_mat(mat_dim)
    μ2, Σ2 = randn(mat_dim), random_psd_mat(mat_dim)

    # Generate samples
    samples1 = rand(MvNormal(μ1, Σ1), n_sample)
    samples2 = rand(MvNormal(μ2, Σ2), n_sample)

    # Objective function
    function f((x_L1, x_D1, x_L2, x_D2))
        -loglikelihood(MvNormal(μ1, decompress_symmetric(x_L1, x_D1)), samples1) - 
            loglikelihood(MvNormal(μ2, decompress_symmetric(x_L2, x_D2)), samples2)
    end
    model = Model(f)

    mat_x0 = random_psd_mat(mat_dim)
    x0 = [mat_x0[Nonconvex.lowertriangind(mat_x0)], diag(mat_x0)]
    x0 = [copy.(x0)..., copy.(x0)...]

    lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
    lbs = [copy.(lbs)..., copy.(lbs)...]

    ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
    ubs = [copy.(ubs)..., copy.(ubs)...]

    addvar!(model, lbs, ubs)

    add_sd_constraint!(model, x -> sd_constraint(x[1:2]))
    add_sd_constraint!(model, x -> sd_constraint(x[3:4]))

    alg = SDPBarrierAlg(sub_alg=IpoptAlg())
    options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
    result = optimize(model, alg, x0, options = options)

    minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
    _Σ1 = sd_constraint(minimizer[1:2])
    _Σ2 = sd_constraint(minimizer[3:4])

    println("result: \n $result")

    println("minimum: \n $minimum")
    println("minimizer: \n $minimizer")
    println("_Σ1: \n $(_Σ1)")
    println("_Σ2: \n $(_Σ2)")

    println("Σ1: \n $Σ1")
    println("Σ2: \n $Σ2")
    println("abs(_Σ1 - Σ1): \n $(abs.(_Σ1 - Σ1))")
    println("abs(_Σ2 - Σ2): \n $(abs.(_Σ2 - Σ2))")
    println("mean(abs(_Σ1 - Σ1)): \n $(mean(abs.(_Σ1 - Σ1)))")
    println("mean(abs(_Σ2 - Σ2)): \n $(mean(abs.(_Σ2 - Σ2)))")

    @test Σ1 ≈ _Σ1 rtol = 0.1
    @test Σ2 ≈ _Σ2 rtol = 0.1
end

@testset "Different algorithm" begin
    # Randomize groundtruth
    μ, Σ = randn(mat_dim), random_psd_mat(mat_dim)

    # Generate samples
    samples = rand(MvNormal(μ, Σ), n_sample)

    # Passing sd_constraint but not using SDPBarrierOptions
    # Objective function
    function f((x_L, x_D))
        -loglikelihood(MvNormal(μ, decompress_symmetric(x_L, x_D)), samples)
    end
    model = Model(f)

    mat_x0 = random_psd_mat(mat_dim)
    x0 = [mat_x0[Nonconvex.lowertriangind(mat_x0)], diag(mat_x0)]
    lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
    ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
    addvar!(model, lbs, ubs)
    add_sd_constraint!(model, sd_constraint)

    result = optimize(model, IpoptAlg(), x0, options=IpoptOptions(max_iter=1, first_order=true))
end
