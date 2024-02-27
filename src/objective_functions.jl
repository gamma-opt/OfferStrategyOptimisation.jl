using JuMP

## -- OBJECTIVE --

function objective!(model::Model, sets::Sets, prices::Prices, costs::Costs; 
    DA_revenue::Bool=true, reserve_revenue::Bool=true, 
    ID_revenue::Bool=true, imbalance_revenue::Bool=true, CCGT_costs::Bool=true, 
    hydro_costs::Bool=true, WVF_piecewise::Bool=true, return_expression::Bool=false)
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W
    π_S = sets.π_S
    π_E = sets.π_E
    π_W = sets.π_W

    DA = DA_revenue ? @expression(model, sum(π_S[s] * prices.DA[t,s] * y[t,s] for t in T, s in S)) : 0

    reserve = reserve_revenue ? @expression(model,sum(π_S[s] * prices.reserve[t,s] * r[t,s] for t in T, s in S )) : 0

    ID = ID_revenue ? @expression(model, sum(π_S[s] * π_E[e] * prices.ID[t, s, e] * z[t,s,e] for t in T, s in S, e in E)) : 0

    imbalance = imbalance_revenue ? @expression(model, sum( π_S[s] * π_E[e] * π_W[w] * (prices.down[t,s,e,w] * delta_excess[t,s,e,w] - prices.up[t,s,e,w] * delta_deficit[t,s,e,w]) for t in T, s in S, e in E, w in W )) : 0


    CCGT = CCGT_costs ? @expression(model, sum( π_S[s] * π_E[e] * ( costs.CCGT_generation * g_CCGT[t,s,e] 
                            + costs.CCGT_startup * u_CCGT_start[t,s,e] 
                            + costs.CCGT_shutdown * u_CCGT_stop[t,s,e])
                            for t in T, s in S, e in E)) : 0

    # Option of piecewise WVF or simple marginal cost in eur/MWh
    if hydro_costs && WVF_piecewise

        println("\n\n\n******************************\n\n WVF piecewise code not tested!!\n\n******************************")

        hydro_startup = @expression(model, sum( π_S[s] * π_E[e] * costs.Hydro_startup * u_hydro_start[t,s,e] for t in T, s in S, e in E))

        hydro_opportunity_cost = @expression(model, - sum( π_S[s] * π_E[e] * (costs.Hydro_WV[k, 1] * l_hydro_k[t, s, e, k] + costs.Hydro_WV[k, 2])
                            for t = T[end], s in S, e in E, k in K))

    elseif hydro_costs && !WVF_piecewise

        hydro_startup = @expression(model, sum( π_S[s] * π_E[e] * costs.Hydro_startup * u_hydro_start[t,s,e] for t in T, s in S, e in E))

        hydro_opportunity_cost = @expression(model, - sum( π_S[s] * π_E[e] * (costs.Hydro_WV[1] * l_hydro[t, s, e] + costs.Hydro_WV[2]) for t = T[end], s in S, e in E))
    else
        hydro_startup = 0
        hydro_opportunity_cost = 0
    end

    if !return_expression
        # create objective
        @objective(model, Max, DA + reserve + ID + imbalance - CCGT - hydro_startup - hydro_opportunity_cost)
    else
        # return only the expression
        @expression(model, DA + reserve + ID + imbalance - CCGT - hydro_startup - hydro_opportunity_cost)
    end
end


## -- OBJECTIVE FOR DEGENERATE SOLUTIONS --
function objective_degenerate_solutions!(model::Model, sets::Sets, delta_excess, delta_deficit)
    
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W

    @objective(model, Min, sum(delta_excess[t,s,e,w] + delta_deficit[t,s,e,w] for t in T, s in S, e in E, w in W))

end
