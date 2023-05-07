using Logging
using JuMP
using DataFrames
using CSV
using Distributions
using Random
using StatsBase
using Gurobi
using PrettyTables
using Dates

@info("Start time    ("*Dates.format(now(), "HH:MM")*")")
##########################################################################################
@info("Initialising parameters")
## -- Time steps --
nT = 24
T = 1:nT
# Time step length: one hour
dt = 1


## -- Scenarios--
dates = ["20220302", "20220309", "20220316", "20220323", "20220330", "20220406", "20220413", 
"20230222", "20230301", "20230308", "20230315"]

nS = length(dates)
nE = 3  # <=
nW = 3
S = 1:nS
E = 1:nE
W = 1:nW

π_S = fill(1/nS, nS)
π_E = fill(1/nE, nE)
π_W = fill(1/nW, nW)


## -- Price steps in offers --
pI = [-500, 0, 25, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 300, 400, 4000.0]
nI = length(pI)
I = 1:nI

pJ = [1.0, 2, 3, 4, 5, 7, 10, 101]
nJ = length(pJ)
J = 1:nJ

@info("S = "*string(nS)*", E = "*string(nE)*", W = "*string(nW)*", I = "*string(nI)*", J = "*string(nJ))


## -- Intraday trading limit (percentage of capacity) --
α = 0.05


## -- Cost and operational parameters --
include("costs_and_operational_parameters.jl")

##########################################################################################

## -- Setting path for results --
path = "RESULTS/I"*string(nI)*"_J"*string(nJ)*"_E"*string(nE)*"/"
isdir(path) || mkdir(path)
@info("path = "*path)

# Results saving
draw_plots = true
save_curves = true
save_results_breakdown = true

##########################################################################################
@info("Generating scenarios")

## -- Scenario generation --
include("Data_DA_and_imbalance/DA_scenario_generation.jl")
include("Data_DA_and_imbalance/imbalance_scenario_generation.jl")
include("Data_FCRD/FCRD_scenario_generation.jl")
include("Data_ID/ID_scenario_generation.jl")
include("Data_wind/wind_scenario_generation.jl")

# Random seed for ID and wind scenario generation
random_seed = 3

DA_prices = DA_price_scenario_generation("Data_DA_and_imbalance/", dates, nT, T, nS, S)
reserve_prices = reserve_price_scenario_generation("Data_FCRD/", dates, nT, T, nS, S)
ID_prices = ID_price_scenario_generation("Data_ID/", dates, nT, T, nS, S, nE, seed = random_seed)

balance_prices_up, balance_prices_down = imbalance_price_scenario_generation("Data_DA_and_imbalance/", 
                    dates, nT, T, nS, S, nE, E, nW; 
                    factor = 2)

wind_factors = wind_scenario_generation("Data_wind/", dates, nT, T, nS, S, nE, nW, 
                                        stdev = 0.074, 
                                        round_to_digits = 8, 
                                        seed = random_seed)


##########################################################################################
@info("Building variables and objective")

## -- MODEL --
model = Model()

## -- Variables --
include("variables.jl")

## -- Objective --
include("objective_functions.jl")
objective = objective!(model, DA_revenue = true,
                              ID_revenue = true,
                              reserve_revenue = true,
                              imbalance_revenue =true,
                              CCGT_costs = true,
                              hydro_costs = true,
                              WVF_piecewise = false)

## -- Constraints --
@info("Building constraints")
include("constraints.jl")

# Market constratins
DA_offer_curve = DA_offer_curve_constraints!(model, x, y, pI, DA_prices) #(TxS)

x_constraints = enforcing_bid_order!(model, x, T, nI) #((I-1)xT)

first_stage_offer_limit = first_stage_offer_limit_constraints!(model, x, v, C_CCGT, C_hydro, C_wind, dt, nI) #(T)

generation_capacity_limit = generation_capacity_constraints!(model, y, r, C_hydro, C_wind, C_CCGT, dt) #(TxS)

reserve_offer_curve = reserve_offer_constraints!(model, v, r, pJ, reserve_prices) #(TxS)

reserve_allocation = reserve_allocation_constraints!(model, r, r_hydro, r_CCGT, T, S, E) #(TxSxE)

second_stage_offer_limit = second_stage_offer_limit_constraints!(model, y, r, z, C_CCGT, C_hydro, C_wind, dt) #(TxSxE)

ID_trade_ub, ID_trade_lb = intraday_trade_constraints!(model, z, α , C_CCGT, C_hydro, C_wind, dt, T, S, E) #2x(TxSxE)

imbalance, excess_ub, deficit_ub = imbalance_settlement_constraints!(model, delta_excess, delta_deficit, y, z, g_CCGT, g_hydro, g_wind) #3x(TxSxExW)

# Wind generation constraints
wind_generation = wind_generation_constraints!(model, g_wind, wind_factors, C_wind, dt, T, S, E, W) #(TxSxExW)

# Hydro generation constraints
hydro_on_off, hydro_start, hydro_stop = hydro_on_off_constraints!(model, u_hydro, u_hydro_start, u_hydro_stop, U_hydro_initial, T, S, E) #3x(TxSxE)

generation_conversion, discharge_lb, discharge_ub = hydro_generation_constraints!(model, g_hydro, f_hydro, u_hydro, r_hydro, conversion_eta, F_min, F_max, dt, T, S, E) #3x(TxSxE)

level_lb, level_ub, level_intertemporal = hydro_water_level_constraints!(model, f_hydro, f_hydro_spill, l_hydro, F_inflow, L_initial, L_min, L_max, N, T, S, E) #3x(TxSxE)

# CCGT operation constraints
CCGT_generation_ub, CCGT_generation_lb = CCGT_generation_limit_constraints!(model, g_CCGT, r_CCGT, u_CCGT, C_CCGT, G_min_level, dt, T, S, E) #2x(TxSxE)

CCGT_on_off, CCGT_start, CCGT_stop = CCGT_on_off_constraints!(model, u_CCGT, u_CCGT_start, u_CCGT_stop, U_CCGT_initial, T, S, E) #3x(TxSxE)

CCGT_min_off_time, CCGT_min_on_time = CCGT_on_off_min_time_constraints!(model, u_CCGT, u_CCGT_start, u_CCGT_stop, U_min_on_time, U_min_off_time, nT, T, S, E) #2x(TxSxE)

ramping_lb, ramping_ub = CCGT_ramping_constraints!(model, g_CCGT, u_CCGT_start, u_CCGT_stop, C_CCGT, G_initial, G_min_level, dt, A_ramp_factor, T, S, E) #3x(TxSxE)

CCGT_reserve_capacity = CCGT_reserve_capacity_constraints!(model, r_CCGT, R_CCGT, T, S, E) #(TxSxE)



##########################################################################################
@info("Solving model    ("*Dates.format(now(), "HH:MM")*")")

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "IntFeasTol"  => 1e-6,
    "TimeLimit"   => 300,
    "LogFile" => path*"solve_1.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)

## -- Solving --
optimize!(model)


first_solve = solve_time(model)
println("solve time: ", first_solve)

##########################################################################################
@info("Optimising over degenerate solutions     ("*Dates.format(now(), "HH:MM")*")")

obj_value = objective_value(model)

obj_expression = objective!(model, DA_revenue = true,
                                ID_revenue = true,
                                reserve_revenue = true,
                                imbalance_revenue =true,
                                CCGT_costs = true,
                                hydro_costs = true,
                                WVF_piecewise = false,
                                return_expression=true)

obj_value_constraint = @constraint(model, obj_expression == obj_value)

new_objective = objective_degenerate_solutions!(model, delta_excess, delta_deficit, T, S, E, W)

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "IntFeasTol"  => 1e-6,
    "MIPFocus"    => 3,
    "TimeLimit"   => 1200,
    "LogFile" => path*"solve_2.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)

## -- Solving --
optimize!(model)



# Saving summary for printing
second_solve = olve_time(model)
println("solve time: ", second_solve)

##########################################################################################
@info("Drawing bidding curves   ("*Dates.format(now(), "HH:MM")*")")

if draw_plots
    include("plot_functions.jl")
    isdir(path*"DA_curves/") || mkdir(path*"DA_curves/")
    isdir(path*"reserve_curves/") || mkdir(path*"reserve_curves/")
    for t in T

        plot_transposed_bid_curve(x[:,t], pI, t, path*"DA_curves/DA"*string(t), 
                    plot_DA_prices=true, default_limits=false, price_limits =[-10, 500],
                    DA_prices=DA_prices[t,:], title = "Hour "*string(t))

        if t == 23
            plot_transposed_reserve_offer_curve(v[:,t], pJ, t, path*"reserve_curves/reserve"*string(t),
                        plot_reserve_prices=true, default_limits=false, price_limits = [-1, 115], 
                        reserve_prices=reserve_prices[t,:], title = "Hour "*string(t))

        else
            plot_transposed_reserve_offer_curve(v[:,t], pJ, t, path*"reserve_curves/reserve"*string(t),
                        plot_reserve_prices=true, default_limits=false, price_limits = [-1, 25], 
                        reserve_prices=reserve_prices[t,:], title = "Hour "*string(t))
        end

    end
end

##########################################################################################
@info("Exporting results")

if save_curves

    ## -- Saving offer curves --
    optimal_DA_offers = DataFrame(Hour = [T...])
    for (i,p) in enumerate(pI)
        col = string(p)
        optimal_DA_offers[!,col] = round.([value.(x[i,:])...], digits = 6)
    end
    CSV.write(path*"optimal_DA_offers.csv", optimal_DA_offers, index=false)

    optimal_reserve_offers = DataFrame(Hour = [T...])
    for (j,p) in enumerate(pJ)
        col = string(p)
        optimal_reserve_offers[!,col] = round.([value.(v[j,:])...], digits = 6)
    end
    CSV.write(path*"optimal_reserve_offers.csv", optimal_reserve_offers, index=false)

end 

##########################################################################################

if save_results_breakdown
    ## -- Saving optimisation results --
    io = open(path*"results.txt", "w");

    write(io, "S = "*string(nS)*"\nE = "*string(nE)*"\nW = "*string(nW)*"\n\nalpha = "*string(α));

    write(io, "\n\n\nI = "*string(nI)*"\nJ = "*string(nJ));

    write(io, "\n\nfirst solution time: ", first_solve)
    write(io, "\nsecond solution time: ", second_solve)

    write(io, "\n\n\n*********************************************\n\n");

    write(io, "Profit Breakdown\n\n")

    profit_DA = sum(π_S[s] * DA_prices[t,s] * value(y[t,s]) for t in T, s in S)
    profit_reserve = sum(π_S[s] * reserve_prices[t,s] * value(r[t,s]) for t in T, s in S)
    profit_ID = sum(π_S[s] * π_E[e] * ID_prices[t,s,e] * value(z[t,s,e]) for t in T, s in S, e in E)
    revenue_ID = sum(π_S[s] * π_E[e] * ID_prices[t,s,e] * (value(z[t,s,e]) >= 0 ? value(z[t,s,e]) : 0)  for t in T, s in S, e in E)
    costs_ID = sum(π_S[s] * π_E[e] * ID_prices[t,s,e] * (value(z[t,s,e]) < 0 ? value(z[t,s,e]) : 0)  for t in T, s in S, e in E)
    profit_imbalances = sum(π_S[s] * π_E[e] * π_W[w] * (balance_prices_down[t,s,e,w] * value(delta_excess[t,s,e,w]) - balance_prices_up[t,s,e,w] * value(delta_deficit[t,s,e,w])) for t in T, s in S, e in E, w in W)
    excess_imbalances_revenue = sum(π_S[s] * π_E[e] * π_W[w] * balance_prices_down[t,s,e,w] * value(delta_excess[t,s,e,w]) for t in T, s in S, e in E, w in W)
    deficit_imbalances_cost = sum(π_S[s] * π_E[e] * π_W[w] * ( - balance_prices_up[t,s,e,w] * value(delta_deficit[t,s,e,w])) for t in T, s in S, e in E, w in W)


    measures =  ["Profit", "DA profit", "Reserve profit", "ID profit", "ID revenue", "ID costs", "Imbalance settlement", "Imbalance revenue", "Imbalance costs"]
    values = [obj_value, profit_DA, profit_reserve, profit_ID, revenue_ID, costs_ID, profit_imbalances, excess_imbalances_revenue, deficit_imbalances_cost]
    profit_breakdown = DataFrame(Measure = measures, Values = values)
    pretty_table(io, profit_breakdown)

    write(io, "\n\n\n*********************************************\n\n\n");

    write(io, "Commitments and operation \n\n")

    commitments = DataFrame(Hour = [T...])
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

    for t in T
        append!(DA_commitment, round(mean(value.(   y[t,:]  )), digits = 4))
        append!(reserve_commitment, round(mean(value.(  r[t,:]  )), digits = 4))
        append!(ID_commitment, round(mean(value.(   z[t,:,:]   )), digits = 4))
        append!(excess, round(mean(value.(  delta_excess[t,:,:,:]   )), digits = 4))
        append!(deficit, round(mean(value.( delta_deficit[t,:,:,:]   )), digits = 4))
        append!(hydro, round(mean(value.(   g_hydro[t,:,:]  )), digits = 4))
        append!(CCGT, round(mean(value.(    g_CCGT[t,:,:]   )), digits = 4))
        append!(wind, round(mean(value.(    g_wind[t,:,:,:] )), digits = 4))
        append!(hydro_reserve, round(mean(value.(   r_hydro[t,:,:]  )), digits = 4))
        append!(CCGT_reserve, round(mean(value.(    r_CCGT[t,:,:]   )), digits = 4))
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
    pretty_table(io, commitments)

    close(io)
end