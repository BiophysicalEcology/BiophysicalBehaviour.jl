# Writes data/reference_thermoregulation.csv. It was run on BiophysicalGeometry 7ff35a4,
# HeatExchange 2085d75 and BiophysicalBehaviour 6970295, with positional shape constructors.

using CSV, DataFrames, BiophysicalGeometry, Unitful
cylinder(mass, b) = Cylinder(mass, 1000.0u"kg/m^3", b)
half_cylinder(mass, b) = HalfCylinder(mass, 1000.0u"kg/m^3", b)
include("reference_cases.jl")

# A case that throws is recorded as one `error` row: there is no result to keep.
rows = mapreduce(vcat, REFERENCE_CASES) do case
    values = try
        reference_values(reference_thermoregulate(case))
    catch err
        @warn "Reference case failed" case err
        ["error" => NaN]
    end
    [(; case..., key, value) for (key, value) in values]
end
CSV.write(joinpath(@__DIR__, "data", "reference_thermoregulation.csv"), DataFrame(rows))
