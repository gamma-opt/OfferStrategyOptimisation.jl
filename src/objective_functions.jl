using JuMP

## -- OBJECTIVE --

function profit_objective(model::Model, sets::Sets, prices::Prices, costs::Costs; 
    DA_revenue::Bool=true, reserve_revenue::Bool=true, 
    ID_revenue::Bool=true, imbalance_revenue::Bool=true, CCGT_costs::Bool=true, 
    hydro_costs::Bool=true, water_value::Bool=true)
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W
    π_S = sets.π_S
    π_E = sets.π_E
    π_W = sets.π_W

    DA = DA_revenue ? @expression(model, sum(π_S[s] * prices.DA[t,s] * model[:y][t,s] for t in T, s in S)) : 0

    reserve = reserve_revenue ? @expression(model,sum(π_S[s] * prices.reserve[t,s] * model[:r][t,s] for t in T, s in S )) : 0

    ID = ID_revenue ? @expression(model, sum(π_S[s] * π_E[e] * prices.ID[t, s, e] * model[:z][t,s,e] for t in T, s in S, e in E)) : 0

    imbalance = imbalance_revenue ? @expression(model, sum( π_S[s] * π_E[e] * π_W[w] * (prices.down[t,s,e,w] * model[:delta_excess][t,s,e,w] - prices.up[t,s,e,w] * model[:delta_deficit][t,s,e,w]) for t in T, s in S, e in E, w in W )) : 0


    CCGT = CCGT_costs ? @expression(model, sum( π_S[s] * π_E[e] * ( costs.CCGT_generation * model[:g_CCGT][t,s,e] 
                            + costs.CCGT_startup * model[:u_CCGT_start][t,s,e] 
                            + costs.CCGT_shutdown * model[:u_CCGT_stop][t,s,e])
                            for t in T, s in S, e in E)) : 0

   

    hydro_startup = hydro_costs ? @expression(model, sum( π_S[s] * π_E[e] * costs.Hydro_startup * model[:u_hydro_start][t,s,e] for t in T, s in S, e in E)) : 0

    hydro_opportunity_cost = water_value ? @expression(model, - sum( π_S[s] * π_E[e] * (costs.Hydro_WV[1] * model[:l_hydro][t, s, e] + costs.Hydro_WV[2]) for t = T[end], s in S, e in E)) : 0


    @expression(model, DA + reserve + ID + imbalance - CCGT - hydro_startup - hydro_opportunity_cost)

end


## -- OBJECTIVE FOR DEGENERATE SOLUTIONS --
function degenerate_solutions_objective(model::Model, sets::Sets)
    
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W

    @expression(model, sum(model[:delta_excess][t,s,e,w] + model[:delta_deficit][t,s,e,w] for t in T, s in S, e in E, w in W))

end
