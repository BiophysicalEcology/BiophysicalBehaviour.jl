module ExampleEnvironments

# The microclimates and the hourly loop of "Get started", for the pages of the manual that need them.

using BiophysicalBehaviour
using Microclimate
using Unitful
using Statistics: mean

export madison_environments, palm_springs_environments, simulate

const CACHE = Ref{Any}(nothing)

"""
    madison_environments()

The available environments of "Get started": the example inputs of Microclimate.jl (the middle day of each month
at Madison, Wisconsin) with no shade and with 90 % shade.
"""
function madison_environments()
    CACHE[] === nothing || return CACHE[]
    model = MicroModel(;
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
    )
    function microclimate(shade)
        inputs = MicroInputs(;
            site = example_site(),
            soil_profile = example_soil_profile(),
            environment_minmax = example_monthly_weather(),
            environment_daily = example_daily_environment(; shade = fill(shade, 12)),
            environment_hourly = example_hourly_environment(),
            initial_soil_moisture = fill(0.1, length(model.depths)),
        )
        solve(MicroProblem(model, inputs))
    end
    CACHE[] = AvailableEnvironments(microclimate(0.0), microclimate(0.9), 0.0, 0.9, model.depths, model.heights)
    return CACHE[]
end

const PALM_SPRINGS = Ref{Any}(nothing)

"""
    palm_springs_environments()

The available environments of the tutorial "A lizard's day": monthly climate at Palm Springs, California, with no
shade and with 90 % shade, and air temperatures at five heights up to 2 m.
"""
function palm_springs_environments()
    PALM_SPRINGS[] === nothing || return PALM_SPRINGS[]
    elevation = 130.0u"m"
    site = example_site(; latitude = 33.8u"°", longitude = -116.5u"°", elevation, albedo = 0.33, roughness_height = 0.0005u"m")
    temperature_min = [3.7, 6.0, 8.2, 11.2, 15.0, 19.1, 22.6, 22.2, 18.8, 13.5, 7.5, 3.5]u"°C"
    temperature_max = [19.3, 21.8, 23.6, 27.0, 31.0, 35.9, 39.2, 38.4, 35.4, 30.3, 23.9, 19.7]u"°C"
    weather = example_monthly_weather(;
        reference_temperature_min = temperature_min,
        reference_temperature_max = temperature_max,
        reference_wind_speed_min = [0.24, 0.26, 0.29, 0.31, 0.31, 0.31, 0.31, 0.29, 0.28, 0.25, 0.24, 0.23]u"m/s",
        reference_wind_speed_max = [2.40, 2.62, 2.91, 3.06, 3.06, 3.06, 3.06, 2.91, 2.76, 2.47, 2.40, 2.33]u"m/s",
        reference_humidity_min = [34.0, 33.7, 33.4, 31.1, 31.4, 31.1, 33.4, 35.2, 34.2, 33.0, 33.1, 34.0] ./ 100,
        reference_humidity_max = [95.6, 94.3, 89.4, 83.3, 82.7, 83.2, 86.2, 89.1, 90.8, 92.2, 94.8, 99.3] ./ 100,
        cloud_cover_min = fill(0.0, 12),
        cloud_cover_max = fill(0.0, 12),
    )
    heights = [1.0, 50.0, 100.0, 150.0, 200.0]u"cm"
    model = MicroModel(; heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
    )
    depths = model.depths
    deep_soil_temperature = mean([mean(u"K".(temperature_min)), mean(u"K".(temperature_max))])
    function microclimate(shade)
        inputs = MicroInputs(;
            site,
            soil_profile = example_soil_profile(depths; bulk_density = 1.5u"Mg/m^3"),
            environment_minmax = weather,
            environment_daily = example_daily_environment(;
                shade = fill(shade, 12),
                rainfall = fill(0.0u"kg/m^2", 12),
                deep_soil_temperature = fill(deep_soil_temperature, 12),
            ),
            environment_hourly = example_hourly_environment(; elevation),
            initial_soil_moisture = fill(0.2, length(depths)),
        )
        solve(MicroProblem(model, inputs))
    end
    PALM_SPRINGS[] = AvailableEnvironments(microclimate(0.0), microclimate(0.9), 0.0, 0.9, depths, heights)
    return PALM_SPRINGS[]
end

"""
    simulate(organism, environments, site, steps)

Call `thermoregulate` for each hour in `steps`, passing the depth chosen in one hour on to the next.
"""
function simulate(organism, environments, site, steps)
    limits = thermoregulation(organism)
    depth = limits.depth.reference
    commenced = false
    map(steps) do step
        (step - 1) % 24 == 0 && (commenced = false)
        out = thermoregulate(organism, environments, limits, site, step, depth; activity_commenced = commenced)
        depth = out.depth_node
        commenced = commenced || !(out.state isa Resting)
        out
    end
end

end
