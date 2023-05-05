using CSV
using DataFrames

function DA_price_scenario_generation(path::String, dates::Vector{String}, nT, T, nS, S)

    if length(dates) == nS

        DA_prices =  zeros(Float64, nT, nS) # dimension (TxS)

        for s in S
            DA_price_df = CSV.read(path*"DA_prices_"*string(dates[s])*".csv", DataFrame)

            for t in T
                hour = filter(row -> row.hour == t, DA_price_df)
                DA_prices[t,s] = hour[1,:price]
            end
        end

        DA_prices
    else
        println("number of scenarios S don't match!!")
    end
end