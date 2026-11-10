using Documenter
using DocumenterFragments: integrate_fragments
using NonconvexCore
using NonconvexIpopt, NonconvexMMA, NonconvexPercival, NonconvexNLopt,
    NonconvexNOMAD, NonconvexMetaheuristics, NonconvexSearch, NonconvexJuniper,
    NonconvexPavito, NonconvexTOBS, NonconvexMultistart, NonconvexSemidefinite,
    NonconvexAugLagLab

# Every package ships its documentation as a DocumenterFragments fragment in
# its `docs/` directory; this site assembles all fragments plus its own pages.
fragdir(pkg) = joinpath(dirname(dirname(pathof(pkg))), "docs")

const FRAGMENTS = [
    (NonconvexCore, "core"),
    (NonconvexIpopt, "ipopt"),
    (NonconvexMMA, "mma"),
    (NonconvexPercival, "percival"),
    (NonconvexNLopt, "nlopt"),
    (NonconvexNOMAD, "nomad"),
    (NonconvexMetaheuristics, "metaheuristics"),
    (NonconvexSearch, "search"),
    (NonconvexJuniper, "juniper"),
    (NonconvexPavito, "pavito"),
    (NonconvexTOBS, "tobs"),
    (NonconvexMultistart, "multistart"),
    (NonconvexSemidefinite, "semidefinite"),
    (NonconvexAugLagLab, "auglag"),
]

main_src = joinpath(@__DIR__, "src")
integration = integrate_fragments(
    main_src,
    [(; dir = fragdir(pkg), mount = slug) for (pkg, slug) in FRAGMENTS],
)

core_fragments = [f for f in integration.fragments if f.name == "NonconvexCore"]
algorithm_fragments = [f for f in integration.fragments if f.name != "NonconvexCore"]

makedocs(
    sitename = "Nonconvex.jl",
    modules = integration.modules,
    pages = Any[
        "Getting started" => "index.md",
        "Models, gradients and results" =>
            Any[[f.pages for f in core_fragments]...],
        "Algorithms" => Any[
            "Overview" => "algorithms/algorithms.md",
            [f.pages for f in algorithm_fragments]...,
            "algorithms/surrogate.md",
        ],
    ],
    plugins = integration.plugins,
    warnonly = true,
    # Fragment code samples are illustrative rather than self-contained, so
    # they are not run as doctests in the assembled site.
    doctest = false,
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/JuliaNonconvex/Nonconvex.jl.git",
        push_preview = true,
    )
end
