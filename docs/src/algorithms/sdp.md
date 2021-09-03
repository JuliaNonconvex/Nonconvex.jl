# Semidifinite programming

## Description

If you have matrix in your optimization objective and you need to keep that positive semidefinite, Nonconvex supports a series of interfaces for [Barrier method for semidefinite programming](http://eaton.math.rpi.edu/faculty/Mitchell/courses/matp6640/notes/24A_SDPbarrierbeamer.pdf), which is implemented as a meta-algorithm transferring your optimization target to a general linear programming problem, then solving by pre-specified `sub_options`.

## Quick start

Optimizating over a multivariate gaussian distribution with artificial samples using `Ipopt`:

```julia
    # Draw random multivariate gaussian samples
    # Random groundtruth
    μ, Σ = randn(mat_dim), random_psd_mat(mat_dim)
    # Generate
    samples = rand(MvNormal(μ, Σ), n_sample)

    # Define objective function
    function f((x_L, x_D))
        -loglikelihood(MvNormal(μ, decompress_symmetric(x_L, x_D)), samples)
    end

    # Define settings
    model = Model(f)
    mat_x0 = random_psd_mat(mat_dim)
    x0 = [mat_x0[Nonconvex.lowertriangind(mat_x0)], diag(mat_x0)]
    lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
    ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
    addvar!(model, lbs, ubs)
    add_sd_constraint!(model, sd_constraint)
    alg = SDPBarrierAlg(sub_alg=IpoptAlg())
    options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))

    # Optimize
    result = optimize(model, alg, x0, options = options)

    # Check result
    minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
    _Σ = sd_constraint(minimizer)
    println("result: \n $result")
    println("minimum: \n $minimum")
    println("minimizer: \n $minimizer")
    println("_Σ: \n $(_Σ)")
    println("Σ: \n $Σ")
    println("abs(_Σ - Σ): \n $(abs.(_Σ - Σ))")
    println("mean(abs(_Σ - Σ)): \n $(mean(abs.(_Σ - Σ)))")

```

## Options

`SDPBarrierOptions` only transfers your optimization target to a regular NLP problem, please specify a `sub_options` to solve that 

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
