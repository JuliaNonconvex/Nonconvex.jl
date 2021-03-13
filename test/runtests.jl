using SafeTestsets

@safetestset "General" begin include("general.jl") end

@safetestset "MMA approximation" begin include("approximation.jl") end

@safetestset "MMA algorithm" begin include("mma.jl") end

@safetestset "AugLag" begin include("auglag.jl") end
