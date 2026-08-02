using BiophysicalBehaviour
using Aqua
using SafeTestsets
using Test

@testset "Quality assurance" begin
    Aqua.test_unbound_args(BiophysicalBehaviour)
    Aqua.test_stale_deps(BiophysicalBehaviour)
    Aqua.test_undefined_exports(BiophysicalBehaviour)
    Aqua.test_project_extras(BiophysicalBehaviour)
    Aqua.test_deps_compat(BiophysicalBehaviour)
end

@safetestset "endotherm" begin include("endotherm.jl") end
@safetestset "ectotherm" begin include("ectotherm.jl") end
@safetestset "transient" begin include("transient.jl") end
@safetestset "trans_behav_r" begin include("trans_behav_r.jl") end
@safetestset "transient_endotherm" begin include("transient_endotherm.jl") end
@safetestset "arrest" begin include("arrest.jl") end
