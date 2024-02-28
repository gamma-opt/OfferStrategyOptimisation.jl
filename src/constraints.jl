using JuMP

## -- CONSTRAINTS --

# -- Day-ahead offer curve clearing (1) --
function DA_offer_curve_constraints!(model::Model, x, y, pI::Vector{Float64}, DA_price::Matrix{Float64})
    offer_curve_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()

    for t in T, s in S
        i = findfirst(pI .≥ DA_price[t,s]) - 1
        c = @constraint(model, y[t,s] == (DA_price[t,s] - pI[i])/(pI[i+1] - pI[i])* x[i+1,t] + (pI[i+1] - DA_price[t,s])/(pI[i+1] - pI[i]) *x[i,t])
        offer_curve_constraints[t,s] = c
    end
    offer_curve_constraints

end 

# Optional constraints for enforcing x_i-1 ≤ x_i
function enforcing_bid_order!(model, x, T, nI)
    x_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for i = 2:nI, t in T
        x_constraints[i, t] = @constraint(model, x[i-1, t] ≤ x[i, t])
    end
    x_constraints
end

# -- Reserve offers clearing (1) --
function reserve_offer_constraints!(model::Model, v, r, pJ::Vector{Float64}, reserve_price::Matrix{Float64})
    offer_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for t in T, s in S
        j_cleared = findfirst(pJ .> reserve_price[t,s]) - 1
        c = @constraint(model, r[t,s] == sum(v[j, t] for j in 1:j_cleared))
        offer_constraints[t,s] = c
    end
    offer_constraints
end

# -- Reserve allocation to hydro and CCGT (1) --
function reserve_allocation_constraints!(model::Model, r, r_hydro, r_CCGT, T, S, E)
    allocation_constraints = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
        allocation_constraints[t,s,e] = @constraint(model, r[t,s] == r_hydro[t,s,e] + r_CCGT[t,s,e])
    end

    allocation_constraints
end

# -- Capacity limits on first stage offers (1) --
function first_stage_offer_limit_constraints!(model::Model, x, v, C_CCGT, C_hydro, C_wind, dt, nI)

    DA_reserve_offer_limit = Dict{Int64, ConstraintRef}()
    for t in T
        c = @constraint(model, x[nI, t] + sum(v[j, t] for j in J)*dt ≤ (C_CCGT + C_hydro + C_wind)*dt )
        DA_reserve_offer_limit[t] = c
    end

    DA_reserve_offer_limit
end

# -- Capacity limits on second stage offers (1) --
function second_stage_offer_limit_constraints!(model::Model, y, r, z, C_CCGT, C_hydro, C_wind, dt)

    ID_offer_limit = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    for t in T, s in S, e in E
        c = @constraint(model, y[t,s] + r[t,s]*dt + z[t,s,e] ≤ (C_CCGT + C_hydro + C_wind)*dt )
        ID_offer_limit[t, s, e] = c
    end
    ID_offer_limit
end

## -- Intraday trade quantity constraint (1) --
function intraday_trade_constraints!(model::Model, z, α , C_CCGT, C_hydro, C_wind, dt, T, S, E)

    trade_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    trade_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
        trade_ub[t,s,e] = @constraint(model, z[t,s,e] ≤ α * (C_CCGT + C_hydro + C_wind) * dt)
        trade_lb[t,s,e] = @constraint(model, z[t,s,e] ≥ - α * (C_CCGT + C_hydro + C_wind) * dt)
    end

    trade_ub, trade_lb
end

# Capacity limit on generation and reserve. Not optional!
function generation_capacity_constraints!(model, y, r, C_hydro, C_wind, C_CCGT, dt)
    
    generation_capacity_limit = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for t in T, s in S

        generation_capacity_limit[t,s] = @constraint(model, y[t,s] + r[t,s] * dt ≤ (C_CCGT + C_wind + C_hydro)*dt)

    end
    generation_capacity_limit
end


# -- Imbalance constraints (3) --
function imbalance_settlement_constraints!(model::Model, delta_excess, delta_deficit, y, z, g_CCGT, g_hydro, g_wind)
    imbalances = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()
    positive_limit = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()
    negative_limit = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E, w in W
        c1 = @constraint(model, delta_excess[t,s,e,w] - delta_deficit[t,s,e,w] == g_CCGT[t,s,e] + g_hydro[t,s,e] + g_wind[t,s,e,w] - (y[t,s] + z[t,s,e]))
        imbalances[t,s,e,w] = c1

        c2 = @constraint(model, delta_excess[t,s,e,w] ≤ (C_CCGT + C_hydro)*dt + g_wind[t,s,e,w])
        positive_limit[t,s,e,w] = c2

        c3 = @constraint(model, delta_deficit[t,s,e,w] ≤ (C_CCGT + C_hydro + C_wind)*dt)
        negative_limit[t,s,e,w] = c3
    end

    imbalances, positive_limit, negative_limit
end


# -- Hydro on-off state constraints (3) --
function hydro_on_off_constraints!(model::Model, u_hydro, u_hydro_start, u_hydro_stop, U_hydro_initial, T, S, E)
    hydro_on_off =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    hydro_start =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    hydro_stop =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
        if t == 1
            hydro_on_off[t,s,e] = @constraint(model, u_hydro[t,s,e] == U_hydro_initial + u_hydro_start[t,s,e] - u_hydro_stop[t,s,e])
        else 
            hydro_on_off[t,s,e] = @constraint(model, u_hydro[t,s,e] == u_hydro[t-1,s,e] + u_hydro_start[t,s,e] - u_hydro_stop[t,s,e])
        end

        hydro_start[t,s,e] = @constraint(model, u_hydro[t,s,e] ≥ u_hydro_start[t,s,e])

        hydro_stop[t,s,e] = @constraint(model, u_hydro[t,s,e] ≤ 1- u_hydro_stop[t,s,e])

    end
    hydro_on_off, hydro_start, hydro_stop
end


# -- Hydro generation constraints (3) --
function hydro_generation_constraints!(model::Model, g_hydro, f_hydro, u_hydro, r_hydro, conversion_eta, F_min, F_max, dt, T, S, E)

    generation_conversion = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    discharge_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    discharge_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()


    for t in T, s in S, e in E
        generation_conversion[t,s,e] = @constraint(model, g_hydro[t,s,e] == conversion_eta * f_hydro[t,s,e] * dt)

        discharge_lb[t,s,e] = @constraint(model, u_hydro[t,s,e] * F_min ≤ f_hydro[t,s,e])#+ r_hydro[t,s,e]/conversion_eta)
        discharge_ub[t,s,e] = @constraint(model, f_hydro[t,s,e] + r_hydro[t,s,e]/conversion_eta ≤ u_hydro[t,s,e] * F_max)
        
    end

    generation_conversion, discharge_lb, discharge_ub

end


# -- Hydro water level constraints (3) --
function hydro_water_level_constraints!(model::Model, f_hydro, f_hydro_spill, l_hydro, F_inflow, L_initial, L_min, L_max, N, T, S, E)

    level_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    level_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    level_intertemporal = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E

        level_lb[t,s,e] = @constraint(model, L_min ≤ l_hydro[t,s,e])
        level_ub[t,s,e] = @constraint(model, l_hydro[t,s,e] ≤ L_max)

        if t == 1
            level_intertemporal[t,s,e] = @constraint(model, l_hydro[t,s,e]- L_initial == N * (- f_hydro[t,s,e] - f_hydro_spill[t,s,e] + F_inflow[t]))
        else
            level_intertemporal[t,s,e] = @constraint(model, l_hydro[t,s,e]- l_hydro[t-1,s,e] == N * (- f_hydro[t,s,e] - f_hydro_spill[t,s,e] + F_inflow[t]))
        end
    end

   level_lb, level_ub, level_intertemporal

end

# -- Conventional generation on-off state constraints (3) --
function CCGT_on_off_constraints!(model::Model, u_CCGT, u_CCGT_start, u_CCGT_stop, U_CCGT_initial, T, S, E)
    CCGT_on_off =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    CCGT_start =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    CCGT_stop =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()


    for t in T, s in S, e in E
        if t == 1
            CCGT_on_off[t,s,e] = @constraint(model, u_CCGT[t,s,e] == U_CCGT_initial + u_CCGT_start[t,s,e] - u_CCGT_stop[t,s,e] )
        else 
            CCGT_on_off[t,s,e] = @constraint(model, u_CCGT[t,s,e] == u_CCGT[t-1,s,e] + u_CCGT_start[t,s,e] - u_CCGT_stop[t,s,e])
        end

        CCGT_start[t,s,e] = @constraint(model, u_CCGT[t,s,e] ≥ u_CCGT_start[t,s,e])

        CCGT_stop[t,s,e] = @constraint(model, u_CCGT[t,s,e] ≤ 1- u_CCGT_stop[t,s,e])
    end

    CCGT_on_off, CCGT_start, CCGT_stop
end

# -- CCGT generation on-off minimum time constraints (3) --
function CCGT_on_off_min_time_constraints!(model::Model, u_CCGT, u_CCGT_start, u_CCGT_stop, U_min_on_time, U_min_off_time, nT, T, S, E)
    
    CCGT_min_off_time = Dict{Tuple{Tuple{Int64, Int64}, Int64, Int64}, ConstraintRef}()
    CCGT_min_on_time = Dict{Tuple{Tuple{Int64, Int64}, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
        for t1 in t:min(t+U_min_off_time, nT)
            CCGT_min_off_time[(t, t1), s, e] = @constraint(model, u_CCGT[t1,s,e] ≤ (1 - u_CCGT_stop[t,s,e]))
        end
        
        for t2 in t:min(t+U_min_on_time, nT)
            CCGT_min_on_time[(t, t2), s, e] = @constraint(model, u_CCGT_start[t,s,e] ≤ u_CCGT[t2,s,e])
        end
    end

    CCGT_min_off_time, CCGT_min_on_time
end

# -- CCGT generation limit constraints (2) --
function CCGT_generation_limit_constraints!(model::Model, g_CCGT, r_CCGT, u_CCGT, C_CCGT, G_min_level, dt, T, S, E)
    generation_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    generation_lb =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
   

    for t in T, s in S, e in E
        generation_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] + r_CCGT[t,s,e] * dt ≤ u_CCGT[t,s,e] * C_CCGT * dt)
        generation_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] ≥ u_CCGT[t,s,e] * C_CCGT * dt * G_min_level)
    end

    generation_ub, generation_lb
end


# -- CCGT generation ramping constraints (2) --
function CCGT_ramping_constraints!(model::Model, g_CCGT, u_CCGT_start, u_CCGT_stop, C_CCGT, G_initial, G_min_level, dt, A_ramp_factor, T, S, E)
    
    ramping_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    ramping_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
    
        if t == 1
            ramping_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] - G_initial ≥ - A_ramp_factor * C_CCGT * dt - u_CCGT_stop[t,s,e] * C_CCGT * dt * (G_min_level- A_ramp_factor))
            ramping_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] - G_initial ≤ A_ramp_factor * C_CCGT * dt + u_CCGT_start[t,s,e] * C_CCGT * dt * (G_min_level- A_ramp_factor))
        else
            ramping_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] - g_CCGT[t-1,s,e] ≥ - A_ramp_factor * C_CCGT * dt - u_CCGT_stop[t,s,e] * C_CCGT * dt * (G_min_level- A_ramp_factor))
            ramping_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] - g_CCGT[t-1,s,e] ≤ A_ramp_factor * C_CCGT * dt + u_CCGT_start[t,s,e] * C_CCGT * dt * (G_min_level- A_ramp_factor))
        end
    end

    ramping_lb, ramping_ub
end


# -- CCGT maximum reserve capability --
function CCGT_reserve_capacity_constraints!(model::Model, r_CCGT, R_CCGT, T, S, E)
    reserve_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E
        reserve_ub[t,s,e] = @constraint(model, r_CCGT[t,s,e] ≤ R_CCGT)
    end
    reserve_ub
end

# -- Wind generation operational constraint (1) --
function wind_generation_constraints!(model::Model, g_wind, W_realisation, C_wind, dt, T, S, E, W)
    wind_generation = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()

    for t in T, s in S, e in E, w in W
        wind_generation[t,s,e,w] = @constraint(model, g_wind[t,s,e,w] ≤ W_realisation[t,s,e,w] * C_wind * dt)
    end
    wind_generation
end

