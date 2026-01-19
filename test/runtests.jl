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
