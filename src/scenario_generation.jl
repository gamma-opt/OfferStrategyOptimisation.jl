using Distributions
using CSV
using DataFrames
using Random
using StatsBase



"Generate DA price scenarios from data files in CSV format including headers hour (in CET), price. The file should be named path*indentifier*date.csv. One scenario is generated per each data point (i.e. date)."
function DA_price_scenario_generation(path::String, indentifier::String, dates::Vector{String}, sets::Sets)
    T = sets.T
    S = sets.S
    nT = sets.T[end]
    nS = sets.S[end]


    if length(dates) == nS

        DA_prices =  zeros(Float64, nT, nS) # dimension (TxS)

        for s in S
            DA_price_df = CSV.read(path*indentifier*string(dates[s])*".csv", DataFrame)

            for t in T
                hour = filter(row -> row.hour == t, DA_price_df)
                DA_prices[t,s] = hour[1,:price]
            end
        end

        DA_prices
    else
        throw(DomainError("Number of scenarios S does not match the dates given!"))
    end
end

"Generate FCRD reserve price scenarios from data files in CSV format. The CSV must include headers hour (in CET) and price. The file should be named path*FCR_up_*date.csv. One scenario is generated per each data point (i.e. date)."
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


"Generate ID reserve price scenarios from data files in CSV format. The CSV must include headers DeliveryHour (in CET) and Price. The file should be named path*ID_Finland_*date.csv. Scenarios are sampled from the available data for that hour during that day."
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



"Generate imbalance settlement price scenarios from data files in CSV format. The CSV must include headers hour (in CET) and difference, which is the different between DA and imbalance price (DA-imbalance) . The file should be named path*imbalance_prices_*date.csv. Scenarios are generated in a rule based manner according to the deviation (i.e. the difference)."
function imbalance_price_scenario_generation(path::String, dates::Vector{String}, nT, T, nS, S, nE, E, nW; factor::Int=2)

    if length(dates) == nS

        price_up = zeros(Float64, nT, nS, nE, nW); # dimension (TxSxExW)     
        price_down = zeros(Float64, nT, nS, nE, nW); # dimension (TxSxExW)

        for s in S
            balance_price_data = CSV.read(path*"imbalance_prices_"*dates[s]*".csv", DataFrame)

            for t in T

                data = filter(row -> row.hour == t, balance_price_data)
                difference = round.(data.difference/factor, digits = 2)
                deviation = abs.(difference)


                # Up-regulation hour. price_down = DA price
                if difference < [0]
                    
                    # up-regulating hour as it is in data
                    low_up = data.price_up - deviation # notice low_up > DA price
                    low_down = data.DA_price

                    mid_up = data.price_up # realisation
                    mid_down = data.DA_price 

                    high_up = data.price_up + deviation # notice low_up >> DA price
                    high_down = data.DA_price
            
    
                # Down-regulation hour. price_up = DA price
                elseif deviation > [0] 
            
                    # down-regulating hour as it is in data
                    low_up = data.DA_price
                    low_down = data.price_down - deviation # notice low_down << DA price

                    mid_up = data.DA_price
                    mid_down = data.price_down # realisation

                    high_up = data.DA_price
                    high_down = data.price_down + deviation # notice high_down < DA price


                # No regulation. price_up = price_down = DA_price
                else 
            
                    low_up = data.DA_price
                    low_down = data.DA_price
                    
                    mid_up = data.DA_price 
                    mid_down = data.DA_price
                    
                    high_up = data.DA_price
                    high_down = data.DA_price
                end

                for e in E
                    price_up[t,s,e,:] = [low_up mid_up high_up]
                    price_down[t,s,e,:] = [low_down mid_down high_down]
                end
            end
        end

        # Retun prices
        price_up, price_down
   
    else
        println("number of scenarios S don't match!!")
    end
end



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