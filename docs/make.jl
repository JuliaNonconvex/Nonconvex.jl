using Documenter, Nonconvex
Nonconvex.@load MMA Ipopt NLopt Percival Bayesian Hyperopt Juniper Pavito MTS Semidefinite

makedocs(
    sitename="Nonconvex.jl",
    pages = [
        "Getting started" => "index.md",
        "Problem definition" => "problem/problem.md",
        "Gradients, Jacobians and Hessians" => "gradients/gradients.md",
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
