using PrettyTables
using Dates
using DataFrames
using CSV
using JLD

function results_breakdown(model::Model, sets::Sets, prices::Prices, costs::Costs, profit_objective::Float64, path::String, run_path::String;
    α = 0.05, first_solve::Float64 = -1.0, first_gap::Float64 = -1.0, second_solve::Float64 = -1.0, second_gap::Float64 = -1.0)


    

    io = open(path*"results_"*run_path*".txt", "w");

    write(io, "Run: "*string(Dates.now())*"\n\n")

    write(io, "S = "*string(sets.S[end])*"\nE = "*string(sets.E[end])*"\nW = "*string(sets.W[end])*"\n\nalpha = "*string(α));

    write(io, "\n\n\nI = "*string(sets.I[end])*"\nJ = "*string(sets.J[end]));
    write(io, "\n\npI = "*string(sets.pI)*"\npJ = "*string(sets.pJ));


    write(io, "\n\n\n*********************************************\n\n");

    write(io, "\n\nfirst solution time: "*string(round(first_solve, digits=2))*"\nfirst gap: "*string(first_gap))
    write(io, "\nsecond solution time: "*string(round(second_solve,digits=2))*"\nsecond gap: "*string(second_gap))

    write(io, "\n\n\n*********************************************\n\n");


    write(io, "Profit Breakdown\n\n")


    profit_DA = sum(sets.π_S[s] * prices.DA[t,s] * value(model[:y][t,s]) for t in sets.T, s in sets.S)

    profit_reserve = sum(sets.π_S[s] * prices.reserve[t,s] * value(model[:r][t,s]) for t in sets.T, s in sets.S)

    profit_ID = sum(sets.π_S[s] * sets.π_E[e] * prices.ID[t,s,e] * value(model[:z][t,s,e]) for t in sets.T, s in sets.S, e in sets.E)

    revenue_ID = sum(sets.π_S[s] * sets.π_E[e] * prices.ID[t,s,e] * ( value(model[:z][t,s,e]) >= 0 ? value(model[:z][t,s,e]) : 0 )  
                    for t in sets.T, s in sets.S, e in sets.E)

    costs_ID = sum(sets.π_S[s] * sets.π_E[e] * prices.ID[t,s,e] * ( value( model[:z][t,s,e]) < 0 ? value(model[:z][t,s,e]) : 0 )  
                    for t in sets.T, s in sets.S, e in sets.E)

    profit_imbalances = sum(sets.π_S[s] * sets.π_E[e] * sets.π_W[w] 
                    * (prices.down[t,s,e,w] * value(model[:delta_excess][t,s,e,w]) - prices.up[t,s,e,w] * value(model[:delta_deficit][t,s,e,w])) 
                    for t in sets.T, s in sets.S, e in sets.E, w in sets.W)

    excess_imbalances_revenue = sum(sets.π_S[s] * sets.π_E[e] * sets.π_W[w] * prices.down[t,s,e,w] * value(model[:delta_excess][t,s,e,w]) 
                    for t in sets.T, s in sets.S, e in sets.E, w in sets.W)

    deficit_imbalances_cost = sum(sets.π_S[s] * sets.π_E[e] * sets.π_W[w] * ( - prices.up[t,s,e,w] * value(model[:delta_deficit][t,s,e,w])) 
                    for t in sets.T, s in sets.S, e in sets.E, w in sets.W)

    water_value = sum( sets.π_S[s] * sets.π_E[e] * (costs.Hydro_WV[1] * value(model[:l_hydro][24, s, e]) + costs.Hydro_WV[2]) 
                    for s in sets.S, e in sets.E)

    operating_costs = - sum(sets.π_S[s] * sets.π_E[e] * costs.CCGT_generation * value(model[:g_CCGT][t,s,e]) for t in sets.T, s in sets.S, e in sets.E)


    measures =  ["Profit (obj value)", "DA profit", "Reserve profit", "ID profit", "ID revenue", "ID costs", "Imbalance settlement", "Imbalance revenue", "Imbalance costs", "Water value", "Operating costs"]

    measure_values = [profit_objective, profit_DA, profit_reserve, profit_ID, revenue_ID, costs_ID, profit_imbalances, excess_imbalances_revenue, deficit_imbalances_cost, water_value, operating_costs]

    profit_breakdown = DataFrame(Measure = measures, Values = measure_values)

    pretty_table(io, profit_breakdown)

    write(io, "\n\n\n*********************************************\n\n\n");

    write(io, "Commitments and operation \n\n")

    commitments = DataFrame(Hour = [sets.T...])
    DA_commitment = []
    reserve_commitment = []
    ID_commitment = []
    excess = []
    deficit = []
    hydro = []
    hydro_reserve = []
    CCGT = []
    CCGT_reserve = []
    wind = []
    wind_curtailment = []

    for t in sets.T
        append!(DA_commitment, mean(value.(   model[:y][t,:]  )))
        append!(reserve_commitment, mean(value.(  model[:r][t,:]  )))
        append!(ID_commitment, mean(value.(   model[:z][t,:,:]   )))
        append!(excess, mean(value.(  model[:delta_excess][t,:,:,:]   )))
        append!(deficit, mean(value.( model[:delta_deficit][t,:,:,:]   )))
        append!(hydro, mean(value.(   model[:g_hydro][t,:,:]  )))
        append!(CCGT, mean(value.(    model[:g_CCGT][t,:,:]   )))
        append!(wind, mean(value.(    model[:g_wind][t,:,:,:] )))
        append!(hydro_reserve, mean(value.(   model[:r_hydro][t,:,:]  )))
        append!(CCGT_reserve, mean(value.(    model[:r_CCGT][t,:,:]   )))
        append!(wind_curtailment, mean( wind_factors[t,:,:,:]*40  - value.( model[:g_wind][t,:,:,:])))
    end
    commitments.DA_commitment = DA_commitment
    commitments.Reserve_commitment = reserve_commitment
    commitments.ID_commitment = ID_commitment
    commitments.Excess = excess
    commitments.Deficit = deficit
    commitments.Hydro = hydro
    commitments.Reserve_hydro = hydro_reserve
    commitments.CCGT = CCGT
    commitments.Reserve_CCGT = CCGT_reserve
    commitments.Wind = wind
    commitments.Wind_curtailment = wind_curtailment

    pretty_table(io, commitments; formatters=ft_printf("%5.2f"))

     # Save values to CSV
    CSV.write(path*"Generation_"*run_path*".csv", commitments, index=false)

    close(io)
end


function save_bid_curves(model::Model, sets::Sets, path::String)

    x = model[:x]
    v = model[:v]

    optimal_DA_offers = DataFrame(Hour = [sets.T...])

    for (i,p) in enumerate(sets.pI)
        col = string(p)
        optimal_DA_offers[!,col] = round.([value.(x[i,:])...], digits = 6)
    end
    CSV.write(path*"optimal_DA_offers.csv", optimal_DA_offers, index=false)

    optimal_reserve_offers = DataFrame(Hour = [sets.T...])

    for (j,p) in enumerate(sets.pJ)
        col = string(p)
        optimal_reserve_offers[!,col] = round.([value.(v[j,:])...], digits = 6)
    end
    CSV.write(path*"optimal_reserve_offers.csv", optimal_reserve_offers, index=false)

end


function write_market_vars_jld(model::Model, sets::Sets, path::String)

    file = jldopen(path*"variable_values.jld","w")

    
    file["x"] = [value(model[:x][i,t]) for i in sets.I, t in sets.T]
    file["v"] = [value(model[:v][j,t]) for j in sets.J, t in sets.T]
    file["r"] = [value(model[:r][t,s]) for t in sets.T, s in sets.S]
    file["z"] = [value(model[:z][t,s,e]) for t in sets.T, s in sets.S, e in sets.E]

    file["delta_excess"] = [value(model[:delta_excess][t,s,e,w]) for t in sets.T, s in sets.S, e in sets.E, w in sets.W]
    file["delta_deficit"] = [value(model[:delta_deficit][t,s,e,w]) for t in sets.T, s in sets.S, e in sets.E, w in sets.W]

    close(file)
end