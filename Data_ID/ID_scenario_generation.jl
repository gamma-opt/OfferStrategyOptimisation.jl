using CSV
using DataFrames
using Random
using StatsBase


function ID_price_scenario_generation(path::String, dates::Vector{String}, nT, T, nS, S, nE; seed=123)
    ID_scenarios = zeros(nT, nS, nE);
    Random.seed!(seed)

    for s in S
        df_FI = CSV.read(path*"ID_Finland_"*dates[s]*".csv", DataFrame);
    
        for t in T
            h_prices = filter(row -> row.DeliveryHour==t, df_FI)[:,:Price]
            
            ID_scenarios[t,s,:] = sample(h_prices, nE, replace=false)
        end
    end

    ID_scenarios;
end