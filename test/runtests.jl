using DEBmicroTrait
using SafeTestsets

@safetestset "metabolism" begin include("test_metabolism.jl") end
@safetestset "supeca" begin include("test_supeca.jl") end
@safetestset "assimilation" begin include("test_assimilation.jl") end
@safetestset "setup" begin include("test_setup.jl") end
@safetestset "supeca" begin include("test_supeca.jl") end
@safetestset "thermostoichwizard" begin include("test_thermostoichwizard.jl") end
@safetestset "turnover" begin include("test_turnover.jl") end
