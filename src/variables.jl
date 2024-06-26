using JuMP

## -- Market variables --  

"Declare variables x, y, v, r, z, Δ⁺, Δ⁻ related to market functions.

If variables for ID or reserve market are not included in the model, the function returns zero matrices for them."
function market_variables!(model::Model, 
                        sets::Sets;
                        include_reserve::Bool = true,
                        include_ID::Bool = true)

    # Readability
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W
    I = sets.I
    J = sets.J


    # DA offer quantities
    @variable(model, x[I, T] ≥ 0)

    # DA dispatch quantities
    @variable(model, y[T, S] ≥ 0)

    if include_reserve
        # Reserve offer quantities
        @variable(model, v[J, T] ≥ 0)

        # Reserve offer quantities
        @variable(model, r[T, S] ≥ 0)
    else
        # Reserve offer quantities
        @variable(model, v[J, T] == 0)

        # Reserve offer quantities
        @variable(model, r[T, S] == 0)
        @warn("Reserve market variables v and r fixed to zero.")
    end

    if include_ID
        # ID order (trade) quantities
        @variable(model, z[T, S, E])
    else
        @variable(model, z[T, S, E] == 0)
        @warn("Intraday market variables z fixed to zero.")
    end

    # Imbalances
    @variable(model, delta_excess[T, S, E, W] ≥ 0)
    @variable(model, delta_deficit[T, S, E, W] ≥ 0)

    return x, y, v, r, z, delta_excess, delta_deficit
end



# -- Hydro power generation --
"Declare hydropower generation related variables g_hydro, u_hydro, u_hydro_start, u_hydro_stop, f_hydro, f_hydro_spill, l_hydro, r_hydro.

If variables for reserve market or on-off functionality are not included in the model, the function returns zero or one matrices for them.
"
function hydropower_variables!(model::Model, 
    sets::Sets;
    include_reserve::Bool = true,
    include_start_stop = true)


    # Readability
    T = sets.T
    S = sets.S
    E = sets.E


    @variable(model, g_hydro[T, S, E] ≥ 0)
    @variable(model, u_hydro[T, S, E], Bin)

    if include_start_stop
        @variable(model, u_hydro_start[T, S, E], Bin)
        @variable(model, u_hydro_stop[T, S, E], Bin)
    else
        u_hydro_start = nothing
        u_hydro_stop = nothing
        @warn("Variables u_hydro_start and u_hydro_stop excluded from model. Make sure hydro_costs turned off in objective.")
    end


    @variable(model, f_hydro[T, S, E] ≥ 0)
    @variable(model, f_hydro_spill[T, S, E] ≥ 0)
    @variable(model, l_hydro[T, S, E] ≥ 0)

    if include_reserve
        @variable(model, r_hydro[T, S, E] ≥ 0)
    else
        @variable(model, r_hydro[T, S, E] == 0)
        @warn("Reserve market variables r_hydro fixed to zero.")
    end

    return g_hydro, u_hydro, u_hydro_start, u_hydro_stop, f_hydro, f_hydro_spill, l_hydro, r_hydro

end



# -- CCGT generation --
"Declare CCGT generation related variables g_CCGT, u_CCGT, u_CCGT_start, u_CCGT_stop, r_CCGT.

If variables for reserve market are not included in the model, the function returns zero matrices for them."
function CCGT_generation_variables!(model::Model, 
    sets::Sets;
    include_reserve::Bool = true)

    # Readability
    T = sets.T
    S = sets.S
    E = sets.E

    @variable(model, g_CCGT[T, S, E] ≥ 0)
    @variable(model, u_CCGT[T, S, E], Bin)
    @variable(model, u_CCGT_start[T, S, E], Bin)
    @variable(model, u_CCGT_stop[T, S, E], Bin)

    if include_reserve
        @variable(model, r_CCGT[T, S, E] ≥ 0)
    else
        @variable(model, r_CCGT[T, S, E] == 0)
        @warn("Reserve market variables r_CCGT fixed to zero.")
    end

    return g_CCGT, u_CCGT, u_CCGT_start, u_CCGT_stop, r_CCGT
end



# -- Wind generation --
"Declare wind generation related variables g_wind."
function wind_generation_variables!(model::Model, sets::Sets)

    # Readability
    T = sets.T
    S = sets.S
    E = sets.E
    W = sets.W
    
    @variable(model, g_wind[T, S, E, W] ≥ 0)

    return g_wind
end