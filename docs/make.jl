using Documenter
using DocumenterFragments: integrate_fragments
using Nonconvex
using NonconvexCore
Nonconvex.@load MMA Ipopt NLopt Percival Bayesian Hyperopt Juniper Pavito MTS Semidefinite Metaheuristics NOMAD TOBS AugLag2

# The package documentation fragments assembled into this site. The docs of
# each package live in that package's repository (docs/ with a fragment.toml);
# the sections covered by the fragments are not duplicated here — these docs
# only retain the sections that no fragment covers (getting started, the
# algorithm overview table and the gradient-modifier pages that predate
# NonconvexCore's docs).
const FRAGMENT_PACKAGES = (
    NonconvexCore,
    NonconvexMMA,
    NonconvexIpopt,
    NonconvexNLopt,
    NonconvexPercival,
    NonconvexBayesian,
    NonconvexMultistart,
    NonconvexJuniper,
    NonconvexPavito,
    NonconvexSearch,
    NonconvexSemidefinite,
    NonconvexMetaheuristics,
    NonconvexNOMAD,
    NonconvexTOBS,
    NonconvexAugLagLab,
)

const MAIN_SRC = joinpath(@__DIR__, "src")
const BUILD = joinpath(@__DIR__, "build")

const MOUNTS = Dict(
    NonconvexCore => "core",
    NonconvexMMA => "algorithms/mma",
    NonconvexIpopt => "algorithms/ipopt",
    NonconvexNLopt => "algorithms/nlopt",
    NonconvexPercival => "algorithms/auglag",
    NonconvexBayesian => "algorithms/bayesian",
    NonconvexMultistart => "algorithms/hyperopt",
    NonconvexJuniper => "algorithms/juniper",
    NonconvexPavito => "algorithms/pavito",
    NonconvexSearch => "algorithms/mts",
    NonconvexSemidefinite => "algorithms/sdp",
    NonconvexMetaheuristics => "algorithms/metaheuristics",
    NonconvexNOMAD => "algorithms/nomad",
    NonconvexTOBS => "algorithms/tobs",
    NonconvexAugLagLab => "algorithms/auglaglab",
)

specs = [(; dir=joinpath(pkgdir(p), "docs"), mount=MOUNTS[p]) for p in FRAGMENT_PACKAGES]
module_map = Dict(string(nameof(p)) => p for p in FRAGMENT_PACKAGES)

# Read-only depot trees would propagate to the copied pages and prevent the
# integration from prepending its namespace blocks; make the sources
# writable first, and clear any stale copies from previous builds.
for p in FRAGMENT_PACKAGES
    run(`chmod -R u+w $(joinpath(pkgdir(p), "docs"))`)
end
for mount in values(MOUNTS)
    rm(joinpath(MAIN_SRC, mount); recursive=true, force=true)
end

integration = integrate_fragments(MAIN_SRC, specs; module_map=module_map)

const CORE = [f for f in integration.fragments if f.name == "NonconvexCore"]
const SOLVERS = [f for f in integration.fragments if f.name != "NonconvexCore"]

makedocs(;
    sitename="Nonconvex.jl",
    modules=integration.modules,
    source=MAIN_SRC,
    build=BUILD,
    doctest=false,
    warnonly=true,
    remotes=nothing,
    plugins=[integration.namespacing],
    format=Documenter.HTML(;
        prettyurls=true, edit_link=nothing, repolink=nothing, inventory_version=""
    ),
    pages=Any[
        "Getting started" => "index.md",
        "Core" => [f.pages for f in CORE],
        "Gradients, Jacobians and Hessians" => [
            "Overview" => "gradients/gradients.md",
            "gradients/user_defined.md",
            "gradients/other_ad.md",
            "gradients/chainrules_fd.md",
            "gradients/history.md",
        ],
        "Algorithms" => Any[
            "Overview" => "algorithms/algorithms.md", [f.pages for f in SOLVERS]...
        ],
    ],
)
