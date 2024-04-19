using JuMP

## -- CONSTRAINTS --

# -- Day-ahead offer curve clearing (1) --
function DA_offer_curve_constraints!(model::Model, sets::Sets, prices::Prices)
    offer_curve_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()

    pI = sets.pI
    x = model[:x]
    y = model[:y]

    for t in sets.T, s in sets.S
        i = findfirst(pI .≥ prices.DA[t,s]) - 1
        c = @constraint(model, y[t,s] == (prices.DA[t,s] - pI[i])/(pI[i+1] - pI[i])* x[i+1,t] + (pI[i+1] - prices.DA[t,s])/(pI[i+1] - pI[i]) *x[i,t])
        offer_curve_constraints[t,s] = c
    end
    offer_curve_constraints

end 

# -- Constraints for enforcing x_i-1 ≤ x_i (1) --
function bid_quantity_order_constraints!(model, sets::Sets)
    x = model[:x]

    x_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for i = 2:sets.I[end], t in sets.T
        x_constraints[i, t] = @constraint(model, x[i-1, t] ≤ x[i, t])
    end
    x_constraints
end

# -- Reserve offers clearing (1) --
function reserve_offer_curve_constraints!(model::Model, sets::Sets, prices::Prices)
    v = model[:v]
    r = model[:r]
    pJ = sets.pJ

    offer_constraints = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for t in sets.T, s in sets.S
        j_cleared = findfirst(pJ .> prices.reserve[t,s]) - 1
        c = @constraint(model, r[t,s] == sum(v[j, t] for j in 1:j_cleared))
        offer_constraints[t,s] = c
    end
    offer_constraints
end

# -- Reserve allocation to hydro and CCGT (1) --
function reserve_allocation_constraints!(model::Model, sets::Sets)
    r = model[:r]
    r_hydro = model[:r_hydro]
    r_CCGT = model[:r_CCGT]

    allocation_constraints = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        allocation_constraints[t,s,e] = @constraint(model, r[t,s] == r_hydro[t,s,e] + r_CCGT[t,s,e])
    end

    allocation_constraints
end

# -- Capacity limits on first stage offers (1) --
function first_stage_capacity_limit_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    x = model[:x]
    v = model[:v]
    J = sets.J
    nI = sets.I[end]

    DA_reserve_offer_limit = Dict{Int64, ConstraintRef}()
    for t in sets.T
        c = @constraint(model, x[nI, t] + sum(v[j, t] for j in J)*params.dt ≤ (params.C_CCGT + params.C_Hydro + params.C_Wind)*params.dt )
        DA_reserve_offer_limit[t] = c
    end

    DA_reserve_offer_limit
end

# -- Capacity limits on second stage offers (1) --
function second_stage_capacity_limit_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    y = model[:y]
    r = model[:r]
    z = model[:z]

    ID_offer_limit = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    for t in sets.T, s in sets.S, e in sets.E
        c = @constraint(model, y[t,s] + r[t,s]*params.dt + z[t,s,e] ≤ (params.C_CCGT + params.C_Hydro + params.C_Wind)*params.dt )
        ID_offer_limit[t, s, e] = c
    end
    ID_offer_limit
end

## -- Intraday trade quantity constraint (1) --
function ID_trade_limit_constraints!(model::Model, sets::Sets, params::OperationalParameters; α = 0.05)
    z = model[:z]

    trade_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    trade_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        trade_ub[t,s,e] = @constraint(model, z[t,s,e] ≤ α * (params.C_CCGT + params.C_Hydro + params.C_Wind) * params.dt)
        trade_lb[t,s,e] = @constraint(model, z[t,s,e] ≥ - α * (params.C_CCGT + params.C_Hydro + params.C_Wind) * params.dt)
    end

    trade_ub, trade_lb
end

# Capacity limit on generation and reserve. Not optional!
function generation_capacity_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    y = model[:y]
    r = model[:r]

    generation_capacity_limit = Dict{Tuple{Int64, Int64}, ConstraintRef}()
    for t in sets.T, s in sets.S

        generation_capacity_limit[t,s] = @constraint(model, y[t,s] + r[t,s] * params.dt ≤ (params.C_CCGT + params.C_Wind + params.C_Hydro) * params.dt)

    end
    generation_capacity_limit
end


# -- Imbalance constraints (3) --
function imbalance_calculation_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    y = model[:y]
    z = model[:z]
    delta_excess = model[:delta_excess]
    delta_deficit = model[:delta_deficit]
    g_CCGT = model[:g_CCGT]
    g_hydro = model[:g_hydro]
    g_wind = model[:g_wind]

    imbalances = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()
    excess_ub = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()
    deficit_ub = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E, w in sets.W
        c1 = @constraint(model, delta_excess[t,s,e,w] - delta_deficit[t,s,e,w] == g_CCGT[t,s,e] + g_hydro[t,s,e] + g_wind[t,s,e,w] - (y[t,s] + z[t,s,e]))
        imbalances[t,s,e,w] = c1

        c2 = @constraint(model, delta_excess[t,s,e,w] ≤ (params.C_CCGT + params.C_Hydro)*params.dt + g_wind[t,s,e,w])
        excess_ub[t,s,e,w] = c2

        c3 = @constraint(model, delta_deficit[t,s,e,w] ≤ (params.C_CCGT + params.C_Hydro + params.C_Wind)*params.dt)
        deficit_ub[t,s,e,w] = c3
    end

    imbalances, excess_ub, deficit_ub
end


# -- Hydro on-off state constraints (3 constraint types) --
function hydro_on_off_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    u_hydro = model[:u_hydro]
    u_hydro_start = model[:u_hydro_start]
    u_hydro_stop = model[:u_hydro_stop]    


    hydro_on_off =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    hydro_start =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    hydro_stop =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        if t == 1
            hydro_on_off[t,s,e] = @constraint(model, u_hydro[t,s,e] == params.U_hydro_initial + u_hydro_start[t,s,e] - u_hydro_stop[t,s,e])
        else 
            hydro_on_off[t,s,e] = @constraint(model, u_hydro[t,s,e] == u_hydro[t-1,s,e] + u_hydro_start[t,s,e] - u_hydro_stop[t,s,e])
        end

        hydro_start[t,s,e] = @constraint(model, u_hydro[t,s,e] ≥ u_hydro_start[t,s,e])

        hydro_stop[t,s,e] = @constraint(model, u_hydro[t,s,e] ≤ 1- u_hydro_stop[t,s,e])

    end
    hydro_on_off, hydro_start, hydro_stop
end


# -- Hydro generation constraints (3 constraint types) --
function hydro_generation_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    g_hydro = model[:g_hydro]
    f_hydro = model[:f_hydro]
    u_hydro = model[:u_hydro]
    r_hydro = model[:r_hydro]

    generation_conversion = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    discharge_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    discharge_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()


    for t in sets.T, s in sets.S, e in sets.E
        generation_conversion[t,s,e] = @constraint(model, g_hydro[t,s,e] == params.eta * f_hydro[t,s,e] * params.dt)

        discharge_lb[t,s,e] = @constraint(model, u_hydro[t,s,e] * params.F_min ≤ f_hydro[t,s,e])
        discharge_ub[t,s,e] = @constraint(model, f_hydro[t,s,e] + r_hydro[t,s,e]/params.eta ≤ u_hydro[t,s,e] * params.F_max)
        
    end

    generation_conversion, discharge_lb, discharge_ub

end


# -- Hydro water level constraints (3) --
function hydro_water_level_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    l_hydro = model[:l_hydro]
    f_hydro = model[:f_hydro]
    f_hydro_spill = model[:f_hydro_spill]

    level_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    level_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    level_intertemporal = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E

        level_lb[t,s,e] = @constraint(model, params.L_min ≤ l_hydro[t,s,e])
        level_ub[t,s,e] = @constraint(model, l_hydro[t,s,e] ≤ params.L_max)

        if t == 1
            level_intertemporal[t,s,e] = @constraint(model, l_hydro[t,s,e]- params.L_initial == params.N * (- f_hydro[t,s,e] - f_hydro_spill[t,s,e] + params.F_inflow[t]))
        else
            level_intertemporal[t,s,e] = @constraint(model, l_hydro[t,s,e]- l_hydro[t-1,s,e] == params.N * (- f_hydro[t,s,e] - f_hydro_spill[t,s,e] + params.F_inflow[t]))
        end
    end

   level_lb, level_ub, level_intertemporal

end


# -- Conventional generation on-off state constraints (3 constraint types) --
function CCGT_on_off_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    u_CCGT = model[:u_CCGT]
    u_CCGT_start = model[:u_CCGT_start]
    u_CCGT_stop = model[:u_CCGT_stop]

    CCGT_on_off =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    CCGT_start =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    CCGT_stop =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        if t == 1
            CCGT_on_off[t,s,e] = @constraint(model, u_CCGT[t,s,e] == params.U_CCGT_initial + u_CCGT_start[t,s,e] - u_CCGT_stop[t,s,e] )
        else 
            CCGT_on_off[t,s,e] = @constraint(model, u_CCGT[t,s,e] == u_CCGT[t-1,s,e] + u_CCGT_start[t,s,e] - u_CCGT_stop[t,s,e])
        end

        CCGT_start[t,s,e] = @constraint(model, u_CCGT[t,s,e] ≥ u_CCGT_start[t,s,e])

        CCGT_stop[t,s,e] = @constraint(model, u_CCGT[t,s,e] ≤ 1- u_CCGT_stop[t,s,e])
    end

    CCGT_on_off, CCGT_start, CCGT_stop
end

# -- CCGT generation minimum on-off time constraints (3 constraint types) --
function CCGT_min_on_off_time_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    u_CCGT = model[:u_CCGT]
    u_CCGT_start = model[:u_CCGT_start]
    u_CCGT_stop = model[:u_CCGT_stop]

    CCGT_min_off_time = Dict{Tuple{Tuple{Int64, Int64}, Int64, Int64}, ConstraintRef}()
    CCGT_min_on_time = Dict{Tuple{Tuple{Int64, Int64}, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        for t1 in t:min(t + params.U_min_off_time, sets.T[end])
            CCGT_min_off_time[(t, t1), s, e] = @constraint(model, u_CCGT[t1,s,e] ≤ (1 - u_CCGT_stop[t,s,e]))
        end
        
        for t2 in t:min(t + params.U_min_on_time, sets.T[end])
            CCGT_min_on_time[(t, t2), s, e] = @constraint(model, u_CCGT_start[t,s,e] ≤ u_CCGT[t2,s,e])
        end
    end

    CCGT_min_off_time, CCGT_min_on_time
end

# -- CCGT generation capacity limit constraints (2 constraint types) --
function CCGT_generation_capacity_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    g_CCGT = model[:g_CCGT]
    r_CCGT = model[:r_CCGT]
    u_CCGT = model[:u_CCGT]

    generation_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    generation_lb =  Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
   

    for t in sets.T, s in sets.S, e in sets.E
        generation_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] + r_CCGT[t,s,e] * params.dt ≤ u_CCGT[t,s,e] * params.C_CCGT * params.dt)
        generation_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] ≥ u_CCGT[t,s,e] * params.C_CCGT * params.dt * params.G_min_level)
    end

    generation_ub, generation_lb
end


# -- CCGT generation ramping constraints (2 constraint types) --
function CCGT_ramping_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    g_CCGT = model[:g_CCGT]
    u_CCGT_stop = model[:u_CCGT_stop]
    u_CCGT_start = model[:u_CCGT_start]

    ramping_lb = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()
    ramping_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
    
        if t == 1
            ramping_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] - params.G_initial ≥ - params.A_ramp_factor * params.C_CCGT * params.dt - u_CCGT_stop[t,s,e] * params.C_CCGT * params.dt * (params.G_min_level- params.A_ramp_factor))
            ramping_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] - params.G_initial ≤ params.A_ramp_factor * params.C_CCGT * params.dt + u_CCGT_start[t,s,e] * params.C_CCGT * params.dt * (params.G_min_level- params.A_ramp_factor))
        else
            ramping_lb[t,s,e] = @constraint(model, g_CCGT[t,s,e] - g_CCGT[t-1,s,e] ≥ - params.A_ramp_factor * params.C_CCGT * params.dt - u_CCGT_stop[t,s,e] * params.C_CCGT * params.dt * (params.G_min_level- params.A_ramp_factor))
            ramping_ub[t,s,e] = @constraint(model, g_CCGT[t,s,e] - g_CCGT[t-1,s,e] ≤ params.A_ramp_factor * params.C_CCGT * params.dt + u_CCGT_start[t,s,e] * params.C_CCGT * params.dt * (params.G_min_level- params.A_ramp_factor))
        end
    end

    ramping_lb, ramping_ub
end


# -- CCGT maximum reserve capability --
function CCGT_reserve_capacity_constraints!(model::Model, sets::Sets, params::OperationalParameters)
    r_CCGT = model[:r_CCGT]

    reserve_ub = Dict{Tuple{Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E
        reserve_ub[t,s,e] = @constraint(model, r_CCGT[t,s,e] ≤ params.R_CCGT)
    end
    reserve_ub
end

# -- Wind generation operational constraint (1 constraint type) --
function wind_generation_constraints!(model::Model, sets::Sets, params::OperationalParameters, wind_factor::Array{Float64})
    g_wind = model[:g_wind]
    
    wind_generation = Dict{Tuple{Int64, Int64, Int64, Int64}, ConstraintRef}()

    for t in sets.T, s in sets.S, e in sets.E, w in sets.W
        wind_generation[t,s,e,w] = @constraint(model, g_wind[t,s,e,w] ≤ wind_factor[t,s,e,w] * params.C_Wind * params.dt)
    end
    wind_generation
end

