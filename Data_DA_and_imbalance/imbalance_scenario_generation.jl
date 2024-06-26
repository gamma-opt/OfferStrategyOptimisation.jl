using CSV
using DataFrames

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




function EV0_imbalance_price_scenario_generation(path::String, dates::Vector{String}, nT, T, nS, S, nE, E, nW; factor::Int=1, default_deviation::Float64=5.0)

    if length(dates) == nS

        price_up = zeros(Float64, nT, nS, nE, nW); # dimension (TxSxExW)     
        price_down = zeros(Float64, nT, nS, nE, nW); # dimension (TxSxExW)

        for s in S
            balance_price_data = CSV.read(path*"imbalance_prices_"*dates[s]*".csv", DataFrame)

            for t in T

                data = filter(row -> row.hour == t, balance_price_data)
                deviation = abs.(round.(data.difference/factor, digits = 2))


                # No regulation hour
                if deviation == [0]
                    
                    low_up = data.DA_price .+ default_deviation
                    low_down = data.DA_price

                    mid_up = data.DA_price
                    mid_down = data.DA_price

                    high_up = data.DA_price
                    high_down = data.DA_price .- default_deviation
            
    
                # Regulation hour
                else
            
                    low_up = data.DA_price + deviation
                    low_down = data.DA_price

                    mid_up = data.DA_price
                    mid_down = data.DA_price

                    high_up = data.DA_price
                    high_down = data.DA_price - deviation


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