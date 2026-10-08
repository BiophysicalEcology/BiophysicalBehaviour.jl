# Thermoregulation reproduces data/reference_thermoregulation.csv, recorded by
# generate_reference.jl. Cases recorded as `error` threw then, so they only need to run now.

using CSV, DataFrames, Test, BiophysicalGeometry, Unitful
cylinder(mass, b) = Cylinder(; mass, density = 1000.0u"kg/m^3", axis_ratio_b = b)
half_cylinder(mass, b) = HalfCylinder(; mass, density = 1000.0u"kg/m^3", axis_ratio_b = b)
include("reference_cases.jl")

reference = CSV.read(joinpath(@__DIR__, "data", "reference_thermoregulation.csv"), DataFrame)
case_columns = collect(string.(keys(first(REFERENCE_CASES))))

for (case, expected) in zip(REFERENCE_CASES, groupby(reference, case_columns; sort = false))
    @testset "$case" begin
        @test all(k -> string(expected[1, k]) == string(case[k]), keys(case))
        out = reference_thermoregulate(case)
        expected.key == ["error"] && continue
        values = Dict(reference_values(out))
        for (key, value) in zip(expected.key, expected.value)
            @test haskey(values, key) && isapprox(values[key], value; rtol = 1e-6, atol = 1e-12)
        end
    end
end
