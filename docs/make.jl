import Hyperopt
using Documenter, Nonconvex

makedocs(
    sitename="Nonconvex.jl",
    pages = [
        "Getting started" => "index.md",
        "Problem definition" => "problem.md",
        "Algorithms" => [
            "Overview" => "algorithms/algorithms.md",
            "algorithms/mma.md",
            "algorithms/ipopt.md",
            "algorithms/nlopt.md",
            "algorithms/auglag.md",
            "algorithms/minlp.md",
            "algorithms/hyperopt.md",
            "algorithms/surrogate.md",
        ],
    ],
)
