using SafeTestsets, Test

@safetestset "Model" begin include("model.jl") end
@testset "MMA" begin
    @safetestset "Approximation" begin include("approximation.jl") end
    @safetestset "Algorithm" begin include("mma.jl") end
end
@safetestset "AugLag/Percival" begin include("percival.jl") end
@safetestset "AugLag2" begin include("auglag.jl") end
@safetestset "Ipopt" begin include("ipopt.jl") end
@safetestset "NLopt" begin include("nlopt.jl") end
@safetestset "Utilities" begin include("utilities.jl") end
@safetestset "Hyperopt" begin include("hyperopt.jl") end
