using Documenter, Nonconvex
Nonconvex.@load MMA Ipopt NLopt Percival Bayesian Hyperopt Juniper Pavito MTS Semidefinite

makedocs(
    sitename="Nonconvex.jl",
    pages = [
        "Getting started" => "index.md",
        "Problem definition" => [
            "Overview" => "problem/problem.md",
            "`Model` definition" => "problem/model.md",
            "`DictModel` definition" => "problem/dict_model.md",
            "Querying models" => "problem/queries.md",
        ],
        "Gradients, Jacobians and Hessians" => [
            "Overview" => "gradients/gradients.md",
            "gradients/user_defined.md",
            "gradients/other_ad.md",
            "gradients/chainrules_fd.md",
            "gradients/sparse.md",
            "gradients/symbolic.md",
            "gradients/implicit.md",
            "gradients/history.md",
        ],
        "Algorithms" => [
            "Overview" => "algorithms/algorithms.md",
            "algorithms/mma.md",
            "algorithms/ipopt.md",
            "algorithms/nlopt.md",
            "algorithms/auglag.md",
            "algorithms/minlp.md",
            "algorithms/hyperopt.md",
            "algorithms/surrogate.md",
            "algorithms/mts.md",
            "algorithms/sdp.md",
            "algorithms/metaheuristics.md",
            "algorithms/nomad.md",
        ],
        "Optimization result" => "result.md"
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/JuliaNonconvex/Nonconvex.jl.git",
        push_preview=true,
    )
end
