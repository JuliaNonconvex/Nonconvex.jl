using Test, SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

using Test, Nonconvex, Pkg

if GROUP == "All" || GROUP == "NonconvexIpopt"
    @safetestset "NonconvexIpopt" begin
        @test_throws ArgumentError using NonconvexIpopt
        Nonconvex.@load Ipopt
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexNLopt"
    @safetestset "NonconvexNLopt" begin
        @test_throws ArgumentError using NonconvexNLopt
        Nonconvex.@load NLopt
        NLoptAlg(:LD_MMA)
    end
elseif GROUP == "All" || GROUP == "NonconvexJuniper"
    @safetestset "NonconvexJuniper" begin
        @test_throws ArgumentError using NonconvexJuniper
        Nonconvex.@load Juniper
        JuniperIpoptAlg()
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMMA"
    @safetestset "NonconvexMMA" begin
        @test_throws ArgumentError using NonconvexMMA
        Nonconvex.@load MMA
        MMA87()
        MMA02()
    end
elseif GROUP == "All" || GROUP == "NonconvexPavito"
    @safetestset "NonconvexPavito" begin
        @test_throws ArgumentError using NonconvexPavito
        Nonconvex.@load Pavito
        PavitoIpoptCbcAlg()
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexBayesian"
    @safetestset "NonconvexBayesian" begin
        @test_throws ArgumentError using NonconvexBayesian
        Nonconvex.@load Bayesian
        BayesOptAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexAugLagLab"
    @safetestset "NonconvexAugLagLab" begin
        @test_throws ArgumentError using NonconvexAugLagLab
        Nonconvex.@load AugLag2
        AugLag2()
    end
elseif GROUP == "All" || GROUP == "NonconvexSemidefinite"
    @safetestset "NonconvexSemidefinite" begin
        @test_throws ArgumentError using NonconvexSemidefinite
        Nonconvex.@load Semidefinite
        SDPBarrierAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexSearch"
    @safetestset "NonconvexSearch" begin
        @test_throws ArgumentError using NonconvexSearch
        Nonconvex.@load Search
        MTSAlg()
        LS1Alg()
    end
elseif GROUP == "All" || GROUP == "NonconvexTOBS"
    @safetestset "NonconvexTOBS" begin
        @test_throws ArgumentError using NonconvexTOBS
        Nonconvex.@load TOBS
        TOBSAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMetaheuristics"
    @safetestset "NonconvexMetaheuristics" begin
        @test_throws ArgumentError using NonconvexMetaheuristics
        Nonconvex.@load Metaheuristics
        MetaheuristicsAlg(ECA)
    end
elseif GROUP == "All" || GROUP == "NonconvexNOMAD"
    @safetestset "NonconvexNOMAD" begin
        @test_throws ArgumentError using NonconvexNOMAD
        Nonconvex.@load NOMAD
        NOMADAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMultistart"
    @safetestset "NonconvexMultistart" begin
        @test_throws ArgumentError using NonconvexMultistart
        Nonconvex.@load Multistart
        HyperoptAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexPercival"
    @safetestset "NonconvexPercival" begin
        @test_throws ArgumentError using NonconvexPercival
        Nonconvex.@load AugLag
        AugLag()
    end
else
    throw(ArgumentError("Unknown test group: $GROUP"))
end
