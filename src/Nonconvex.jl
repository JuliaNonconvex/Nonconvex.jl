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
        AugLag,
        IpoptAlg,
        NLoptAlg,
        PercivalAlg,
        KKTCriteria,
        IpoptCriteria,
        FunctionWrapper,
        MMAOptions,
        AugLagOptions,
        IpoptOptions,
        NLoptOptions,
        PercivalOptions,
        Tolerance

using Parameters, Zygote, ChainRulesCore, ForwardDiff
using Ipopt, NLopt, Percival
using LinearAlgebra, Setfield, Requires, SparseArrays, Reexport
using Optim: Optim, AbstractOptimizer
@reexport using LinearAlgebra

abstract type Workspace end

# General

include("utilities/params.jl")
include("functions/functions.jl")
include("functions/value_jacobian.jl")
include("functions/function_docs.jl")

# Models

include("models/model.jl")
include("models/model_docs.jl")
include("mma_approximation/mma_approx.jl")
include("mma_approximation/xmma_approx.jl")
include("mma_approximation/mma_approx_docs.jl")
include("models/mma_model.jl")
include("models/dual_model.jl")

# MMA

include("utilities/convergence.jl")
include("utilities/callbacks.jl")
include("utilities/options.jl")
include("algorithms/mma_algorithm.jl")

# MMA-AugLag

include("functions/aggregations.jl")
include("algorithms/stoch_optimizers.jl")
include("algorithms/nonstoch_optimizers.jl")

# Augmented Lagrangian

include("models/auglag_model.jl")
include("algorithms/auglag_algorithm.jl")

# Ipopt

include("wrappers/ipopt.jl")
include("wrappers/nlopt.jl")
include("wrappers/percival.jl")

end
