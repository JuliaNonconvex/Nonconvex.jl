module Nonconvex

const debugging = Ref(false)
const show_residuals = Ref(false)

export  Model,
        DictModel,
        addvar!,
        add_ineq_constraint!,
        add_eq_constraint!,
        getmin,
        getmax,
        isinteger,
        setmin!,
        setmax!,
        setinteger!,
        optimize,
        Workspace,
        MMA87,
        MMA02,
        MMALag,
        IpoptAlg,
        NLoptAlg,
        AugLag,
        PercivalAlg,
        JuniperIpoptAlg,
        KKTCriteria,
        IpoptCriteria,
        FunctionWrapper,
        MMAOptions,
        IpoptOptions,
        NLoptOptions,
        AugLagOptions,
        PercivalOptions,
        JuniperIpoptOptions,
        Tolerance,
        @constructor

using Parameters, Zygote, ChainRulesCore, ForwardDiff
using Ipopt, NLopt, ADNLPModels, Percival, NLPModelsModifiers
import MathOptInterface, JuMP
const MOI = MathOptInterface
using LinearAlgebra, Setfield, Requires, SparseArrays, Reexport
using Juniper, NamedTupleTools
using Optim: Optim, AbstractOptimizer
using SparseArrays, OrderedCollections
using JuMP: VariableRef, is_binary, is_integer, has_lower_bound,
            has_upper_bound, lower_bound, upper_bound,
            start_value, ConstraintRef, constraint_object,
            AffExpr, objective_function, objective_sense
            
@reexport using LinearAlgebra, OrderedCollections

abstract type Workspace end

# General

include("utilities/params.jl")
include("functions/functions.jl")
include("functions/value_jacobian.jl")
include("functions/function_docs.jl")
include("functions/counting_function.jl")
include("functions/aggregations.jl")

# Models

include("models/flatten.jl")
include("models/model.jl")
include("models/vec_model.jl")
include("models/dict_model.jl")
include("models/model_docs.jl")
include("mma_approximation/mma_approx.jl")
include("mma_approximation/xmma_approx.jl")
include("mma_approximation/mma_approx_docs.jl")
include("models/mma_model.jl")
include("models/dual_model.jl")
include("models/mmalag_model.jl")
include("models/jump.jl")

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

# Wrappers

include("wrappers/ipopt.jl")
include("wrappers/nlopt.jl")
include("wrappers/percival.jl")
include("wrappers/moi.jl")
include("wrappers/juniper.jl")

end
