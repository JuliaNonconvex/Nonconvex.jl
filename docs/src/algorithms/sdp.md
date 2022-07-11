# `NonconvexSemidefinite.jl`

## Description

If you need to keep your any matrix-valued function of the decision variables positive semidefinite, Nonconvex supports an interface for the [barrier method for semidefinite programming](http://eaton.math.rpi.edu/faculty/Mitchell/courses/matp6640/notes/24A_SDPbarrierbeamer.pdf), which is a meta-algorithm transforming the optimization problem to a series of nonlinear programming problems and solving them using the pre-specified `sub_alg` and `sub_options`.

## Quick start

Optimizing over a multivariate gaussian distribution with artificial samples using `Ipopt`:

```julia
using Nonconvex, Distributions
Nonconvex.@load Semidefinite Ipopt

# Draw random multivariate gaussian samples
# Random groundtruth
mat_dim = 3
μ = randn(mat_dim)
Σ = rand(mat_dim, mat_dim)
Σ = Σ + Σ' + 2I
# Generate
n_sample = 1000
samples = rand(MvNormal(μ, Σ), n_sample)

# Define objective function
function f((x_L, x_D))
    return -loglikelihood(MvNormal(μ, decompress_symmetric(x_L, x_D)), samples)
end
# Define the matrix-valued function
function sd_constraint((x_L, x_D))
    return decompress_symmetric(x_L, x_D)
end

# Define settings
model = Model(f)
mat_x0 = rand(mat_dim, mat_dim)
mat_x0 = mat_x0 + mat_x0' + I

x0 = [mat_x0[NonconvexSemidefinite.lowertriangind(mat_x0)], diag(mat_x0)]
lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
addvar!(model, lbs, ubs)
add_sd_constraint!(model, sd_constraint)
alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200), n_iter = 20)

# Optimize
result = optimize(model, alg, x0, options = options)

# Relative error norm
norm(sd_constraint(result.minimizer) - Σ) / norm(Σ)
```

## Options

```@docs
SDPBarrierOptions
```

## Optimizer

```@docs
SDPBarrierAlg
```

## Matrix interfaces

For every `n*n` real positive semidefinite matrix that optimization objective contains, please have two inputs `x_L` and `x_D` representing the lower-triangular and the diagonal part of it. In the function, call `decompress_symmetric(x_L, x_d)` to represent that matrix, which will be handled by Nonocnvex automatically

```@docs
decompress_symmetric
```
