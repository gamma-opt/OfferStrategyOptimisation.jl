using DataFrames
using CSV
using Distributions
using Random

function wind_scenario_generation(path::String,  dates::Vector{String}, nT, T, nS, S, nE, nW; stdev=0.1, round_to_digits=4, seed=123)

    W_scenarios = zeros(nT, nS, nE, nW);
    Random.seed!(seed)
    
    for s in S
    
        df = CSV.read(path*"wind_data_"*dates[s]*".csv", DataFrame)

        for t in T

            # forecast for hour t in scenario s
            realisation = filter(row -> row.hour == t, df)[1, :wind_factor_forecast];

            # Truncated normal distrition [0,1] with mean = wind factor forecast, sd = sample forecast error sd
            dist = truncated(Normal(realisation, stdev), 0.0, 1.0);

            # Sampling ExW scenarios dfrom dist
            W_scenarios[t, s, :, :] = round.(rand(dist, nE, nW), digits = round_to_digits);

        end
    end

    W_scenarios

end