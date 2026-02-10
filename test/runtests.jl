using Test, Nonconvex, Pkg

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "NonconvexIpopt"
    @testset "NonconvexIpopt" begin
        @test_throws ArgumentError using NonconvexIpopt
        Nonconvex.@load Ipopt
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexNLopt"
    @testset "NonconvexNLopt" begin
        @test_throws ArgumentError using NonconvexNLopt
        Nonconvex.@load NLopt
        NLoptAlg(:LD_MMA)
    end
elseif GROUP == "All" || GROUP == "NonconvexJuniper"
    @testset "NonconvexJuniper" begin
        @test_throws ArgumentError using NonconvexJuniper
        Nonconvex.@load Juniper
        JuniperIpoptAlg()
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMMA"
    @testset "NonconvexMMA" begin
        @test_throws ArgumentError using NonconvexMMA
        Nonconvex.@load MMA
        MMA87()
        MMA02()
    end
elseif GROUP == "All" || GROUP == "NonconvexPavito"
    @testset "NonconvexPavito" begin
        @test_throws ArgumentError using NonconvexPavito
        Nonconvex.@load Pavito
        PavitoIpoptCbcAlg()
        IpoptAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexBayesian"
    @testset "NonconvexBayesian" begin
        @test_throws ArgumentError using NonconvexBayesian
        Nonconvex.@load Bayesian
        BayesOptAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexAugLagLab"
    @testset "NonconvexAugLagLab" begin
        @test_throws ArgumentError using NonconvexAugLagLab
        Nonconvex.@load AugLag2
        AugLag2()
    end
elseif GROUP == "All" || GROUP == "NonconvexSemidefinite"
    @testset "NonconvexSemidefinite" begin
        @test_throws ArgumentError using NonconvexSemidefinite
        Nonconvex.@load Semidefinite
        SDPBarrierAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexSearch"
    @testset "NonconvexSearch" begin
        @test_throws ArgumentError using NonconvexSearch
        Nonconvex.@load Search
        MTSAlg()
        LS1Alg()
    end
elseif GROUP == "All" || GROUP == "NonconvexTOBS"
    @testset "NonconvexTOBS" begin
        @test_throws ArgumentError using NonconvexTOBS
        Nonconvex.@load TOBS
        TOBSAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMetaheuristics"
    @testset "NonconvexMetaheuristics" begin
        @test_throws ArgumentError using NonconvexMetaheuristics
        Nonconvex.@load Metaheuristics
        MetaheuristicsAlg(ECA)
    end
elseif GROUP == "All" || GROUP == "NonconvexNOMAD"
    @testset "NonconvexNOMAD" begin
        @test_throws ArgumentError using NonconvexNOMAD
        Nonconvex.@load NOMAD
        NOMADAlg()
    end
elseif GROUP == "All" || GROUP == "NonconvexMultistart"
    @testset "NonconvexMultistart" begin
        @test_throws ArgumentError using NonconvexMultistart
        Nonconvex.@load Multistart
        HyperoptAlg(IpoptAlg())
    end
elseif GROUP == "All" || GROUP == "NonconvexPercival"
    @testset "NonconvexPercival" begin
        @test_throws ArgumentError using NonconvexPercival
        Nonconvex.@load AugLag
        AugLag()
    end
else
    throw(ArgumentError("Unknown test group: $GROUP"))
end
