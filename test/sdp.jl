# Test cases of semidefinite programming
using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences, Distributions
using ChainRulesTestUtils

const FDM = FiniteDifferences

# Test setting
mat_dim = 3
mat_length = mat_dim^2
n_sample = 3000

test_rrule(Nonconvex.rearrange_x, rand((mat_dim^2-mat_dim)÷2), rand(mat_dim))

function random_psd_mat(mat_dim)
    _mat = randn(mat_dim, mat_dim)
    return _mat' * _mat
end

# Randomize groundtruth
μ, Σ = randn(mat_dim), random_psd_mat(mat_dim)
mn = MvNormal(μ, Σ)

# Generate samples
samples = rand(mn, n_sample)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())

options = SDPBarrierOptions(c_init=1, c_decr=0.6, n_iter=5, mat_dim=mat_dim, sub_options=IpoptOptions(max_iter=100, first_order=true))

# Objective function
function f((μ_, x_L, x_D))
    -loglikelihood(MvNormal(μ_, decompress_symmetric(x_L, x_D)), samples)
end

function sd_function((μ_, x_L, x_D))
    decompress_symmetric(x_L, x_D)
end

model = Model(f)

mat_x0 = random_psd_mat(mat_dim)
x0 = (zeros(mat_dim), mat_x0[Nonconvex.lowertriangind(mat_x0)], diag(mat_x0))

addvar!(model, [-Inf*ones(length(_x0)) for _x0 in x0], [Inf*ones(length(_x0)) for _x0 in x0])

add_sd_constraint!(model, MatrixFunctionWrapper(sd_function, mat_dim))

result = optimize(model, alg, x0, options = options)

minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
_μ, _Σ = minimizer[1], sd_function(minimizer)

println("result: \n $result")

println("minimum: \n $minimum")
println("minimizer: \n $minimizer")
println("_μ: \n $_μ")
println("_Σ: \n $(_Σ)")

println("Σ: \n $Σ")
println("μ: \n $μ")
println("abs(_Σ - Σ): \n $(abs.(_Σ - Σ))")
println("abs(_μ - μ): \n $(abs.(_μ - μ))")
println("mean(abs(_Σ - Σ)): \n $(mean(abs.(_Σ - Σ)))")
println("mean(abs(_μ - μ)): \n $(mean(abs.(_μ - μ)))")

@test mean(abs.(_Σ - Σ)) < 0.1 && mean(abs.(_μ - μ)) < 0.1
