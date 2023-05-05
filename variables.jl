using JuMP

## -- VARIABLES --  

# DA offer quantities
@variable(model, x[I, T] ≥ 0)

# DA dispatch quantities
@variable(model, y[T, S] ≥ 0)

# Reserve offer quantities
@variable(model, v[J, T] ≥ 0)

# Reserve offer quantities
@variable(model, r[T, S] ≥ 0)

# ID order (trade) quantities
@variable(model, z[T, S, E])

# Imbalances
@variable(model, delta_excess[T, S, E, W] ≥ 0)
@variable(model, delta_deficit[T, S, E, W] ≥ 0)


# -- Hydro power generation --
@variable(model, g_hydro[T, S, E] ≥ 0)
@variable(model, u_hydro[T, S, E], Bin)
@variable(model, u_hydro_start[T, S, E], Bin)
@variable(model, u_hydro_stop[T, S, E], Bin)
@variable(model, f_hydro[T, S, E] ≥ 0)
@variable(model, f_hydro_spill[T, S, E] ≥ 0)
@variable(model, l_hydro[T, S, E] ≥ 0)
@variable(model, r_hydro[T, S, E] ≥ 0)

# -- Conventional generation --
@variable(model, g_CCGT[T, S, E] ≥ 0)
@variable(model, u_CCGT[T, S, E], Bin)
@variable(model, u_CCGT_start[T, S, E], Bin)
@variable(model, u_CCGT_stop[T, S, E], Bin)
@variable(model, r_CCGT[T, S, E] ≥ 0)

# -- Wind generation --
@variable(model, g_wind[T, S, E, W] ≥ 0)