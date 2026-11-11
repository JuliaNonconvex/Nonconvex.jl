using Test, Nonconvex

@testset "Macro regression tests" begin
    # Unsupported algorithms throw a proper ErrorException (wrapped in a
    # LoadError at macro expansion), not a bare String.
    err = try
        @eval Nonconvex.@load NotAnAlgorithm
        nothing
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ErrorException

    # The generated install_and_load_module code must not leak the modname
    # binding into the caller's scope.
    Nonconvex.install_and_load_module(:NonconvexCore) |> eval
    @test !isdefined(@__MODULE__, :modname)
end

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "NonconvexIpopt"
    @testset "NonconvexIpopt" begin
        @test_throws ArgumentError using NonconvexIpopt
        Nonconvex.@load Ipopt
        IpoptAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexNLopt"
    @testset "NonconvexNLopt" begin
        @test_throws ArgumentError using NonconvexNLopt
        Nonconvex.@load NLopt
        NLoptAlg(:LD_MMA)
    end
end

if GROUP == "All" || GROUP == "NonconvexJuniper"
    @testset "NonconvexJuniper" begin
        @test_throws ArgumentError using NonconvexJuniper
        Nonconvex.@load Ipopt
        Nonconvex.@load Juniper
        JuniperIpoptAlg()
        IpoptAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexMMA"
    @testset "NonconvexMMA" begin
        @test_throws ArgumentError using NonconvexMMA
        Nonconvex.@load MMA
        MMA87()
        MMA02()
    end
end

if GROUP == "All" || GROUP == "NonconvexPavito"
    @testset "NonconvexPavito" begin
        @test_throws ArgumentError using NonconvexPavito
        Nonconvex.@load Ipopt
        Nonconvex.@load Pavito
        PavitoIpoptCbcAlg()
        IpoptAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexBayesian"
    @testset "NonconvexBayesian" begin
        @test_throws ArgumentError using NonconvexBayesian
        Nonconvex.@load Ipopt
        Nonconvex.@load Bayesian
        BayesOptAlg(IpoptAlg())
    end
end

if GROUP == "All" || GROUP == "NonconvexAugLagLab"
    @testset "NonconvexAugLagLab" begin
        @test_throws ArgumentError using NonconvexAugLagLab
        Nonconvex.@load AugLag2
        AugLagAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexSemidefinite"
    @testset "NonconvexSemidefinite" begin
        @test_throws ArgumentError using NonconvexSemidefinite
        Nonconvex.@load Ipopt
        Nonconvex.@load Semidefinite
        SDPBarrierAlg(IpoptAlg())
    end
end

if GROUP == "All" || GROUP == "NonconvexSearch"
    @testset "NonconvexSearch" begin
        @test_throws ArgumentError using NonconvexSearch
        Nonconvex.@load Search
        MTSAlg()
        LS1Alg()
    end
end

if GROUP == "All" || GROUP == "NonconvexTOBS"
    @testset "NonconvexTOBS" begin
        @test_throws ArgumentError using NonconvexTOBS
        Nonconvex.@load TOBS
        TOBSAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexMetaheuristics"
    @testset "NonconvexMetaheuristics" begin
        @test_throws ArgumentError using NonconvexMetaheuristics
        Nonconvex.@load Metaheuristics
        MetaheuristicsAlg(ECA)
    end
end

if GROUP == "All" || GROUP == "NonconvexNOMAD"
    @testset "NonconvexNOMAD" begin
        @test_throws ArgumentError using NonconvexNOMAD
        Nonconvex.@load NOMAD
        NOMADAlg()
    end
end

if GROUP == "All" || GROUP == "NonconvexMultistart"
    @testset "NonconvexMultistart" begin
        @test_throws ArgumentError using NonconvexMultistart
        Nonconvex.@load Ipopt
        Nonconvex.@load Multistart
        HyperoptAlg(IpoptAlg())
    end
end

if GROUP == "All" || GROUP == "NonconvexPercival"
    @testset "NonconvexPercival" begin
        @test_throws ArgumentError using NonconvexPercival
        Nonconvex.@load AugLag
        AugLag()
    end
end

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(Nonconvex)
end

@testset "JET" begin
    using JET
    JET.test_package(Nonconvex; target_modules=(Nonconvex,), analyze_from_definitions=true)
end
