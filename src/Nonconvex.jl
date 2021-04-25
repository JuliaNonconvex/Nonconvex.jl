module Nonconvex

const debugging = Ref(false)
const show_residuals = Ref(false)

export  Model,
        addvar!,
        add_ineq_constraint!,
        add_eq_constraint!,
        getmin,
        getmax,
        setmin!,
        setmax!,
        optimize,
        Workspace,
        MMA87,
        MMA02,
        MMALag,
        IpoptAlg,
        NLoptAlg,
        AugLag,
        PercivalAlg,
        MsOAlg,
        KKTCriteria,
        IpoptCriteria,
        FunctionWrapper,
        MMAOptions,
        IpoptOptions,
        NLoptOptions,
        AugLagOptions,
        MsOOptions,
        PercivalOptions,
        Tolerance,
        get_starting_point,
        MsOResult

using Parameters, Zygote, ChainRulesCore, ForwardDiff
using Ipopt, NLopt, ADNLPModels, Percival, NLPModelsModifiers, MultistartOptimization
using LinearAlgebra, Setfield, Requires, SparseArrays, Reexport
using Optim: Optim, AbstractOptimizer

@reexport using LinearAlgebra

abstract type Workspace end

# General

include("utilities/params.jl")
include("functions/functions.jl")
include("functions/value_jacobian.jl")
include("functions/function_docs.jl")
include("functions/counting_function.jl")
include("functions/aggregations.jl")

# Models

include("models/model.jl")
include("models/model_docs.jl")
include("mma_approximation/mma_approx.jl")
include("mma_approximation/xmma_approx.jl")
include("mma_approximation/mma_approx_docs.jl")
include("models/mma_model.jl")
include("models/dual_model.jl")
include("models/mmalag_model.jl")

# MMA

include("utilities/convergence.jl")
include("utilities/callbacks.jl")
include("utilities/options.jl")
include("algorithms/mma_algorithm.jl")

# MMA-AugLag2

include("algorithms/stoch_optimizers.jl")
include("algorithms/nonstoch_optimizers.jl")
include("algorithms/ammal.jl")

# Augmented Lagrangian

include("models/auglag_model.jl")
include("algorithms/auglag_algorithm.jl")

# Ipopt

include("wrappers/ipopt.jl")
include("wrappers/nlopt.jl")
include("wrappers/percival.jl")
include("wrappers/multistartoptimization.jl")

end
