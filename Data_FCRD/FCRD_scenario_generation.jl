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

#############
# dates = ["20210303", "20210310", "20210317", "20210324", "20210331", 
#         "20220302", "20220309", "20220316", "20220323", "20220330",
#         "20230301", "20230308", "20230315"] 
#reserve_price_scenario_generation("", dates, 24, 1:24, 13, 1:13)