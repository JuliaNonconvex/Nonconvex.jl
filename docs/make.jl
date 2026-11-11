using Documenter
using DocumenterFragments: integrate_fragments
using Nonconvex
using NonconvexCore
Nonconvex.@load MMA Ipopt NLopt Percival Bayesian Hyperopt Juniper Pavito MTS Semidefinite Metaheuristics NOMAD TOBS AugLag2

# The package documentation fragments assembled into this site. The docs of
# each package live in that package's repository (docs/ with a fragment.toml);
# the sections covered by the fragments are not duplicated here — these docs
# only retain the landing page, the algorithm overview table and the @load
# guide that are specific to the Nonconvex.jl umbrella package itself.
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

# Fragment pages cross-link into the core mount with a `../core/...` path
# written for a flat layout in which every fragment page sits directly under
# its mount. The subdirectory mounts used here put each fragment page one
# level deeper, so add the missing parent segment to keep the cross-fragment
# links valid.
for mount in values(MOUNTS)
    for (root, _, files) in walkdir(joinpath(MAIN_SRC, mount))
        for f in files
            endswith(f, ".md") || continue
            path = joinpath(root, f)
            text = read(path, String)
            occursin("../core/", text) || continue
            write(path, replace(text, "../core/" => "../../core/"))
        end
    end
end

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
        prettyurls=true,
        edit_link=nothing,
        repolink=nothing,
        inventory_version="",
        assets=["assets/custom.css"],
    ),
    pages=Any[
        "Home" => "index.md",
        "General" => Any[
            "Loading algorithm packages" => "general/load.md", [f.pages for f in CORE]...
        ],
        "Algorithms" => Any[
            "Overview" => "algorithms/algorithms.md", [f.pages for f in SOLVERS]...
        ],
    ],
)
