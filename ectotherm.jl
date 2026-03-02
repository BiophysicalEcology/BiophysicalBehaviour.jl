using HeatExchange
using Microclimate
using ModelParameters
using Unitful, UnitfulMoles
using Roots
using Test
using Plots
using CSV, DataFrames

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

Tb_NMR = (DataFrame(CSV.File("$testdir/data/TC.csv")))[:, 2] .* u"°C"
# define the geometry
mass = 0.04u"kg"
ρ_body = 1000.0u"kg/m^3"
shapeb = 3
shapec = 2 / 3
# define body shape as an ellipsoid struct of type 'Shape' and give it required values
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) 
# construct a Body, which is naked - this constructor will apply the 'geometry' function 
# to the inputs and return a struct that has the struct for the 'Shape' type, as well 
# as the insulation and the geometry struct
geometric_traits = Body(shape_body, Naked()) 
# construct the Model which holds the parameters of the organism in the Organism concrete struct, 
# of type AbstractOrganism
lizard = Model(Organism(geometric_traits, MorphoPars(), PhysioPars()))

# specify place and time
latitude = -30.0u"°"
longitude = 140.0u"°"
elevation = 10.0u"m"
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [0.01]u"m"
α_substrate = 0.8
albedos = [0.2, 0.2]

# set the environmental parameters
environmental_params = EnvironmentalPars(
    elevation = elevation,
    α_substrate = Param(α_sub, bounds=(0.0, 1.0)),
)

# define daily weather and soil moisture
minima_times = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
maxima_times = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
air_temperature_min = [10.0, 8.0]u"°C" # minimum air temperatures (°C)
air_temperature_max = [30.0, 25]u"°C" # maximum air temperatures (°C)
humidity_min = [20.0, 30.0] # min relative humidity (%)
humidity_max = [80.0, 90.0] # max relative humidity (%)
wind_min = [0.1, 0.2]u"m/s" # min wind speed (m/s)
wind_max = [1.0, 1.4]u"m/s" # max wind speed (m/s)
cloud_min = [20.0, 23.0] # min cloud cover (%)
cloud_max = [90.0, 100.0] # max cloud cover (%)
initial_soil_moisture = [0.2, 0.2] # fractional

n_hours = length(air_temperature_min) * 24

minimum_shade = 0.0
maximum_shade = 0.9

# run microclimate model in minshade environment
micro_min_shade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    minima_times,
    maxima_times,
    initial_soil_moisture,
    albedos,
    shades = fill(minimum_shade, n_hours),
);
plot(1:1:n, u"°C".(micro_min_shade.soil_temperature[:, 1]))

# run microclimate model in minshade environment
micro_max_shade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    minima_times,
    maxima_times,
    initial_soil_moisture,
    albedos,
    shades = fill(maximum_shade, n_hours),
);

min_shade_habitat = EnvironmentalVarsVec(
    air_temperature = K.(micro_min_shade.air_temperature), # second column is first node above surface
    sky_temperature = micro_min_shade.sky_temperature,
    substrate_temperature = micro_min_shade.soil_temperature, # surface temperature
    relative_humidity = Matrix(micro_min_shade.relative_humidity), # second column is first node above surface
    wind_speed = Matrix(micro_min_shade.wind_speed), # second column is first node above surface
    atmospheric_pressure = fill(101325.0u"Pa", n_hours),
    substrate_conductivity = micro_min_shade.soil_thermal_conductivity,
    solar_flux = micro_min_shade.global_solar .* (1.0 - minimum_shade),
    direct_solar_flux = micro_min_shade.direct_solar .* (1.0 - minimum_shade),
    diffuse_solar_flux = micro_min_shade.diffuse_solar .* (1.0 - minimum_shade),
    zenith_angle = micro_min_shade.zenith_angle
);

max_shade_habitat = EnvironmentalVarsVec(
    air_temperature = K.(micro_max_shade.air_temperature), # second column is first node above surface
    sky_temperature = micro_max_shade.sky_temperature,
    substrate_temperature = micro_max_shade.soil_temperature, # surface temperature
    relative_humidity = Matrix(micro_max_shade.relative_humidity), # second column is first node above surface
    wind_speed = Matrix(micro_max_shade.wind_speed), # second column is first node above surface
    atmospheric_pressure = fill(101325.0u"Pa", n_hours),
    substrate_conductivity = micro_max_shade.soil_thermal_conductivity,
    solar_flux = micro_max_shade.global_solar .* (1.0 - maximum_shade),
    direct_solar_flux = micro_max_shade.direct_solar .* (1.0 - maximum_shade),
    diffuse_solar_flux = micro_max_shade.diffuse_solar .* (1.0 - maximum_shade),
    zenith_angle = micro_max_shade.zenith_angle
);

# compute body temperature
shade = 0.0
depth = 1
height = 2
n = n_hours
deep_rh = 99.
deep_vel = 0.01u"m/s"

balances = map(1:n) do i
    atmospheric_pressure = min_shade_habitat.atmospheric_pressure[i]
    zenith_angle     = min_shade_habitat.zenith_angle[i]
    if depth > 1
        air_temperature   = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        sky_temperature   = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        substrate_temperature   = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        relative_humidity      = deep_rh
        wind_speed     = deep_vel
        substrate_conductivity   = min_shade_habitat.substrate_conductivity[i, depth] * (1 - shade) + max_shade_habitat.substrate_conductivity[i, depth] * shade
        solar_flux   = 0.0u"W/m^2"
        direct_solar_flux   = 0.0u"W/m^2"
        diffuse_solar_flux   = 0.0u"W/m^2"
    else
        air_temperature   = min_shade_habitat.air_temperature[i, height + 1] * (1 - shade) + max_shade_habitat.air_temperature[i, height + 1] * shade
        sky_temperature   = min_shade_habitat.sky_temperature[i] * (1 - shade) + max_shade_habitat.sky_temperature[i] * shade
        substrate_temperature   = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        relative_humidity      = min_shade_habitat.relative_humidity[i, height + 1] * (1 - shade) + max_shade_habitat.relative_humidity[i, height + 1] * shade
        wind_speed     = min_shade_habitat.wind_speed[i, height + 1] * (1 - shade) + max_shade_habitat.wind_speed[i, height + 1] * shade
        substrate_conductivity   = min_shade_habitat.substrate_conductivity[i, depth] * (1 - shade) + max_shade_habitat.substrate_conductivity[i, depth] * shade
        solar_flux   = min_shade_habitat.solar_flux[i] * (1 - shade) + max_shade_habitat.solar_flux[i] * shade
        direct_solar_flux   = min_shade_habitat.direct_solar_flux[i] * (1 - shade) + max_shade_habitat.direct_solar_flux[i] * shade
        diffuse_solar_flux   = min_shade_habitat.diffuse_solar_flux[i] * (1 - shade) + max_shade_habitat.diffuse_solar_flux[i] * shade
    end
    env_i = EnvironmentalVars(;
        air_temperature,
        sky_temperature,
        substrate_temperature,
        relative_humidity,
        wind_speed,
        atmospheric_pressure,
        zenith_angle,
        substrate_conductivity,
        solar_flux,
        direct_solar_flux,
        diffuse_solar_flux,
    )
    variables_i = (organism = OrganismalVars(), environment = env_i)
    get_Tb(lizard, environmental_params, variables_i);
end
balance_out = flip2vectors(balances); # pull out each output as a vector
resp_out = flip2vectors(balance_out.resp_out); # pull out each output as a vector
evap_out = flip2vectors(balance_out.evap_out); # pull out each output as a vector
conv_out = flip2vectors(balance_out.conv_out); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.core_temperature), ylims=[-10.0, 55.0])
plot!(1:1:n, micro_min_shade.air_temperature[:, 2])
plot!(1:1:n, micro_min_shade.soil_temperature[:, 1])
plot!(1:1:n, °C.(micro_min_shade.sky_temperature[:, 1]))

function set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
    atmospheric_pressure = min_shade_habitat.atmospheric_pressure[i]
    zenith_angle     = min_shade_habitat.zenith_angle[i]
    if depth > 1
        air_temperature = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        sky_temperature = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        substrate_temperature = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        relative_humidity = deep_rh
        wind_speed = deep_vel
        substrate_conductivity = min_shade_habitat.substrate_conductivity[i, depth] * (1 - shade) + max_shade_habitat.substrate_conductivity[i, depth] * shade
        solar_flux = 0.0u"W/m^2"
        direct_solar_flux = 0.0u"W/m^2"
        diffuse_solar_flux = 0.0u"W/m^2"
    else
        air_temperature = min_shade_habitat.air_temperature[i, height+1] * (1 - shade) + max_shade_habitat.air_temperature[i, height+1] * shade
        sky_temperature = min_shade_habitat.sky_temperature[i] * (1 - shade) + max_shade_habitat.sky_temperature[i] * shade
        substrate_temperature = min_shade_habitat.substrate_temperature[i, depth] * (1 - shade) + max_shade_habitat.substrate_temperature[i, depth] * shade
        relative_humidity = min_shade_habitat.relative_humidity[i, height+1] * (1 - shade) + max_shade_habitat.relative_humidity[i, height+1] * shade
        wind_speed = min_shade_habitat.wind_speed[i, height+1] * (1 - shade) + max_shade_habitat.wind_speed[i, height+1] * shade
        substrate_conductivity = min_shade_habitat.substrate_conductivity[i, depth] * (1 - shade) + max_shade_habitat.substrate_conductivity[i, depth] * shade
        solar_flux = min_shade_habitat.solar_flux[i] * (1 - shade) + max_shade_habitat.solar_flux[i] * shade
        direct_solar_flux = min_shade_habitat.direct_solar_flux[i] * (1 - shade) + max_shade_habitat.direct_solar_flux[i] * shade
        diffuse_solar_flux = min_shade_habitat.diffuse_solar_flux[i] * (1 - shade) + max_shade_habitat.diffuse_solar_flux[i] * shade
    end
    env_i = EnvironmentalVars(;
        air_temperature,
        sky_temperature,
        substrate_temperature,
        relative_humidity,
        wind_speed,
        atmospheric_pressure,
        zenith_angle,
        substrate_conductivity,
        solar_flux,
        direct_solar_flux,
        diffuse_solar_flux,
    )
    return env_i
end

forage_temperature_min = u"K"(24.5u"°C")
forage_temperature_max = u"K"(34.5u"°C")
lethal_temperature_min = u"K"(0.1u"°C")
shade = 0.0
Δshade = 0.01
min_depth = 2
max_height = 2
depth = 1
height = 1
climber = false
burrower = false
shadeseeker = false

Tb = nothing
balances = map(1:n) do i
    #i=12
    height = 1 # reset height
    shade = minimum_shade # reset shade
    depth = 1 # reset depth to surface
    env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
    variables_i = (organism=OrganismalVars(), environment=env_i)
    Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
    # check if too hot and, if so, seek shade
    if shadeseeker
    while Tb > forage_temperature_max && shade <= maximum_shade && depth == 1
        shade += Δshade
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
    end
    shade = clamp(shade, minimum_shade, maximum_shade)
    end
    if climber
    while Tb > forage_temperature_max && (shade >= maximum_shade || !shadeseeker) && height < max_height
        shade = minimum_shade # reset to minimum shade burrow
        height = height + 1
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
        while (Tb > forage_temperature_max || Tb < lethal_temperature_min) && height < max_height
            height = height + 1
            env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
            variables_i = (organism=OrganismalVars(), environment=env_i)
            Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
        end
    end
    end
    # check if too cold or too hot and shade maxed out, and, if so, go underground
    if burrower
    while (Tb < forage_temperature_min || (Tb > forage_temperature_max && shade >= maximum_shade) || (Tb > forage_temperature_max && height >= max_height)) && depth == 1
        shade = minimum_shade # reset to minimum shade burrow
        depth = max(min_depth, depth + 1)
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
        while (Tb > forage_temperature_max || Tb < lethal_temperature_min) && depth < 10
            depth = max(min_depth, depth + 1)
            env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
            variables_i = (organism=OrganismalVars(), environment=env_i)
            Tb = get_Tb(lizard, environmental_params, variables_i).core_temperature;
        end
    end
    end
    return get_Tb(lizard, environmental_params, variables_i);
end
balance_out = flip2vectors(balances); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.core_temperature), ylims=[0, 55.0])
plot!(1:1:n, Tb_NMR;
        xlabel="time", ylabel="body temperature", lw=2,
        linestyle=:dash, linecolor="grey"
    )


################################################

using HeatExchange
using Microclimate
using ModelParameters
using Unitful, UnitfulMoles
using Roots
using Test
using Plots

# define the geometry
mass = 0.04kg
ρ_body = 1000.0kg/m^3
shapeb = 3
shapec = 2 / 3
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) # define body shape as a Cylinder struct of type 'Shape' and give it required values
geometric_traits = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(geometric_traits, MorphoPars(), PhysioPars()))

# specify place and time
latitude = -30.0°
longitude = 140.0°
elevation = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [1.0,]u"cm"
α_substrate = 0.8

# set the environmental parameters
environmental_params = EnvironmentalPars(
    elevation = elevation,
    α_substrate = Param(α_substrate, bounds=(0.0, 1.0)),
)

# define daily weather and soil moisture
minima_times = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
maxima_times = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
air_temperature_min = [10.0, 8.0]u"°C" # minimum air temperatures (°C)
air_temperature_max = [30.0, 25]u"°C" # maximum air temperatures (°C)
humidity_min = [20.0, 30.0] # min relative humidity (%)
humidity_max = [80.0, 90.0] # max relative humidity (%)
wind_min = [0.1, 0.2]u"m/s" # min wind speed (m/s)
wind_max = [1.0, 1.4]u"m/s" # max wind speed (m/s)
cloud_min = [20.0, 23.0] # min cloud cover (%)
cloud_max = [90.0, 100.0] # max cloud cover (%)
initial_soil_moisture = [0.2, 0.2] # fractional
min_shade = 0.0
max_shade = 90.0

# run microclimate model in minshade environment
micro_minshade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    minima_times,
    maxima_times,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    initial_soil_moisture,
    shades = fill(min_shade, length(days)),
)

# run microclimate model in maxshade environment
micro_maxshade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    minima_times,
    maxima_times,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    initial_soil_moisture,
    shades = fill(max_shade, length(days)),
)

env_minshade = EnvironmentalVarsVec(
    air_temperature = K.(micro_minshade.air_temperature[:, 2]), # second column is first node above surface
    sky_temperature = micro_minshade.sky_temperature,
    substrate_temperature = micro_minshade.soil_temperature[:, 1], # surface temperature
    relative_humidity = micro_minshade.relative_humidity[:, 2], # second column is first node above surface
    wind_speed = micro_minshade.wind_speed[:, 2], # second column is first node above surface
    solar_flux = micro_minshade.global_solar .* (1.0 - min_shade / 100.0),
    direct_solar_flux = micro_minshade.direct_solar .* (1.0 - min_shade / 100.0),
    diffuse_solar_flux = micro_minshade.diffuse_solar .* (1.0 - min_shade / 100.0),
    zenith_angle = micro_minshade.zenith_angle
)

env_maxshade = EnvironmentalVarsVec(
    air_temperature = K.(micro_minshade.air_temperature[:, 2]), # second column is first node above surface
    sky_temperature = micro_maxshade.sky_temperature,
    substrate_temperature = micro_maxshade.soil_temperature[:, 1], # surface temperature
    relative_humidity = micro_maxshade.relative_humidity[:, 2], # second column is first node above surface
    wind_speed = micro_maxshade.wind_speed[:, 2], # second column is first node above surface
    solar_flux = micro_maxshade.global_solar .* (1.0 - max_shade / 100.0),
    direct_solar_flux = micro_maxshade.direct_solar .* (1.0 - max_shade / 100.0),
    diffuse_solar_flux = micro_maxshade.diffuse_solar .* (1.0 - max_shade / 100.0),
    zenith_angle = micro_maxshade.zenith_angle
)

# set shade
environment = env_minshade

# compute body temperature
n = length(days) * (length(hours) - 1)
balances = map(1:n) do i
    env_i = EnvironmentalVars(
        air_temperature = environment.air_temperature[i],
        sky_temperature = environment.sky_temperature[i],
        substrate_temperature = environment.substrate_temperature[i],
        relative_humidity = environment.relative_humidity[i],
        wind_speed = environment.wind_speed[i],
        atmospheric_pressure = environment.atmospheric_pressure[i],
        zenith_angle = environment.zenith_angle[i],
        substrate_conductivity = environment.substrate_conductivity[i],
        solar_flux = environment.solar_flux[i],
        direct_solar_flux = environment.direct_solar_flux[i],
        diffuse_solar_flux = environment.diffuse_solar_flux[i],
    )
    variables_i = (organism = OrganismalVars(), environment = env_i)
    get_Tb(lizard, environmental_params, variables_i)
end
balance_out = flip2vectors(balances); # pull out each output as a vector
resp_out = flip2vectors(balance_out.resp_out); # pull out each output as a vector
evap_out = flip2vectors(balance_out.evap_out); # pull out each output as a vector
conv_out = flip2vectors(balance_out.conv_out); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.core_temperature), ylims=[0.0, 55.0])
plot!(1:1:n, environment.air_temperature)

plot(1:1:n, u"mg/hr".(evap_out.m_cut))
plot!(1:1:n, u"mg/hr".(evap_out.m_resp))
plot!(1:1:n, u"mg/hr".(evap_out.m_eyes))
plot!(1:1:n, u"mg/hr".(evap_out.m_evap))