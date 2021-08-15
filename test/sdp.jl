# Test cases of semidefinite programming
using Nonconvex, LinearAlgebra, Test, Zygote, FiniteDifferences, Distributions
const FDM = FiniteDifferences

# Test setting
mat_dim = 5
mat_length = mat_dim^2
n_sample = 300
f_dim = mat_dim+mat_dim*mat_dim

function random_possemi_mat(mat_dim)
    _mat = randn(mat_dim, mat_dim)
    return _mat' * _mat
end

# Randomize groundtruth
real_μ, real_Σ = randn(mat_dim), random_possemi_mat(mat_dim)
real_normal = MvNormal(real_μ, real_Σ)

# Generate samples
samples = rand(real_normal, n_sample)

# Objective function
f = x -> loglikelihood(MvNormal(x[mat_length+1:end], reshape(x[begin:mat_length], (mat_dim, mat_dim))), samples)
obj = MatrixFunctionWrapper(f, mat_dim, f_dim)

# Optimization setting
m = Model(obj)
addvar!(m, [-Inf for _ in 1:f_dim], [Inf for _ in 1:f_dim])
options = SDPOptions(mat_dim, sub_options=IpoptOptions(first_order=true))
alg = SDPAlg(sub_alg=IpoptAlg())

# Add sdp constraints
set_A0(options, 0)
for i in 1:mat_length
    Ai = zeros(mat_dim, mat_dim)
    Ai[CartesianIndices(Ai)[i]] = 1
    add_Ai(options, Ai)
end

x0 = randn(mat_dim^2+mat_dim)
r = Nonconvex.optimize(m, alg, x0, options = options)

print(real_μ)
print(real_Σ)
print(r.minimum)
print(r.minimizer)
print(r.minimizer_mat)
