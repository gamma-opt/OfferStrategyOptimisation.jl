using CSV
using DataFrames

function reserve_price_scenario_generation(path::String, dates::Vector{String}, nT, T, nS, S)

    if length(dates) == nS

        reserve_prices =  zeros(Float64, nT, nS) # dimension (TxS)

        for s in S
            reserve_price_df = CSV.read(path*"FCR_up_"*string(dates[s])*".csv", DataFrame)

            for t in T
                hour = filter(row -> row.hour == t, reserve_price_df)
                reserve_prices[t,s] = hour[1,:price]
            end
        end

        reserve_prices
    else
        println("number of scenarios S don't match!!")
    end
end
