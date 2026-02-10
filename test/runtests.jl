using Test, SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

using Test, Nonconvex, Pkg

if GROUP == "All" || GROUP == "NonconvexIpopt"
    @test_throws ArgumentError using NonconvexIpopt
    Nonconvex.@load Ipopt
    IpoptAlg()
elseif GROUP == "All" || GROUP == "NonconvexNLopt"
    @test_throws ArgumentError using NonconvexNLopt
    Nonconvex.@load NLopt
    NLoptAlg(:LD_MMA)
elseif GROUP == "All" || GROUP == "NonconvexJuniper"
    @test_throws ArgumentError using NonconvexJuniper
    Nonconvex.@load Juniper
    JuniperIpoptAlg()
    IpoptAlg()
elseif GROUP == "All" || GROUP == "NonconvexMMA"
    @test_throws ArgumentError using NonconvexMMA
    Nonconvex.@load MMA
    MMA87()
    MMA02()
elseif GROUP == "All" || GROUP == "NonconvexPavito"
    @test_throws ArgumentError using NonconvexPavito
    Nonconvex.@load Pavito
    PavitoIpoptCbcAlg()
    IpoptAlg()
elseif GROUP == "All" || GROUP == "NonconvexBayesian"
    @test_throws ArgumentError using NonconvexBayesian
    Nonconvex.@load Bayesian
    BayesOptAlg(IpoptAlg())
elseif GROUP == "All" || GROUP == "NonconvexAugLagLab"
    @test_throws ArgumentError using NonconvexAugLagLab
    Nonconvex.@load AugLag2
    AugLag2()
elseif GROUP == "All" || GROUP == "NonconvexSemidefinite"
    @test_throws ArgumentError using NonconvexSemidefinite
    Nonconvex.@load Semidefinite
    SDPBarrierAlg(IpoptAlg())
elseif GROUP == "All" || GROUP == "NonconvexSearch"
    @test_throws ArgumentError using NonconvexSearch
    Nonconvex.@load Search
    MTSAlg()
    LS1Alg()
elseif GROUP == "All" || GROUP == "NonconvexTOBS"
    @test_throws ArgumentError using NonconvexTOBS
    Nonconvex.@load TOBS
    TOBSAlg()
elseif GROUP == "All" || GROUP == "NonconvexMetaheuristics"
    @test_throws ArgumentError using NonconvexMetaheuristics
    Nonconvex.@load Metaheuristics
    MetaheuristicsAlg(ECA)
elseif GROUP == "All" || GROUP == "NonconvexNOMAD"
    @test_throws ArgumentError using NonconvexNOMAD
    Nonconvex.@load NOMAD
    NOMADAlg()
elseif GROUP == "All" || GROUP == "NonconvexMultistart"
    @test_throws ArgumentError using NonconvexMultistart
    Nonconvex.@load Multistart
    HyperoptAlg(IpoptAlg())
elseif GROUP == "All" || GROUP == "NonconvexPercival"
    @test_throws ArgumentError using NonconvexPercival
    Nonconvex.@load AugLag
    AugLag()
else
    throw(ArgumentError("Unknown test group: $GROUP"))
end
