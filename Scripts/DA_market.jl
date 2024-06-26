using JuMP
using Gurobi
using Dates

## -- Setting path for results --
folder = "Results/"
run_path = "DA_market_"*Dates.format(now(), "yyyymmdd-HHMM")
path = folder*run_path*"/"
isdir(path) || mkdir(path)
@info("path = "*path)


include("../src/structures.jl")
include("../src/scenario_generation.jl")
include("../src/variables.jl")
include("../src/objective_functions.jl")
include("../src/constraints.jl")
include("../Plotting/results_analysis.jl")
include("../Plotting/plot_functions.jl")


@info("Start time    ("*Dates.format(now(), "HH:MM")*")")


## -- Data, sets, costs and and operational parameters --

@info("Initialising parameters")

dates = ["20220301", "20220302", "20220308", "20220309", "20220315", "20220316", "20230221","20230222", "20230228", "20230301"]

sets = Sets(S = 1:length(dates), πS = fill(1/length(dates), length(dates)))
costs = Costs()
params = OperationalParameters()


## -- Scenario generation --

@info("Generating scenarios")

prices = Prices()

random_seed = 3

prices.DA = DA_price_scenario_generation("Data_DA_and_imbalance/", "DA_prices_", dates, sets);

prices.reserve = reserve_price_scenario_generation("Data_FCRD/", "FCR_up_", dates, sets);

prices.ID = ID_price_scenario_generation("Data_ID/", "ID_Finland_",  dates, sets, seed = random_seed);

prices.up, prices.down = imbalance_price_scenario_generation("Data_DA_and_imbalance/", "imbalance_prices_", dates, sets; factor = 2);

wind_factors = wind_scenario_generation("Data_wind/", "wind_data_", dates, sets, stdev = 0.0667, round_to_digits = 8, seed = random_seed);


## -- MODEL and variables --

@info("Adding variables and objective to the model")

model = Model()

market_variables!(model, sets, include_reserve = false, include_ID = false);

hydropower_variables!(model, sets, include_start_stop = false, include_reserve = false); 

CCGT_generation_variables!(model, sets, include_reserve = false);

wind_generation_variables!(model, sets);


## -- Objective --

objective = profit_objective!(model, sets, prices, costs,
                        DA_revenue = true,
                        ID_revenue = false,
                        reserve_revenue = false,
                        imbalance_revenue = true,
                        CCGT_costs = true,
                        hydro_costs = false,
                        water_value = true)

@objective(model, Max, objective)


## -- Constraints --

@info("Building constraints")

# Market constratins
DA_offer_curve = DA_offer_curve_constraints!(model, sets, prices)

x_constraints = bid_quantity_order_constraints!(model, sets)


capacity_limit1 = first_stage_capacity_limit_constraints!(model, sets, params)

capacity_limit2 = second_stage_capacity_limit_constraints!(model, sets, params)

α = 0.00

generation_capacity_limit = generation_capacity_constraints!(model, sets, params) 

imbalance, excess_ub, deficit_ub = imbalance_calculation_constraints!(model, sets, params)

# Hydro generation constraints
generation_conversion, discharge_lb, discharge_ub = hydro_generation_constraints!(model, sets, params) 

level_lb, level_ub, level_intertemporal = hydro_water_level_constraints!(model, sets, params)

# CCGT operation constraints
CCGT_on_off, CCGT_start, CCGT_stop = CCGT_on_off_constraints!(model, sets, params)

CCGT_min_off_time, CCGT_min_on_time = CCGT_min_on_off_time_constraints!(model, sets, params)

CCGT_generation_ub, CCGT_generation_lb = CCGT_generation_capacity_constraints!(model, sets, params)

ramping_lb, ramping_ub = CCGT_ramping_constraints!(model, sets, params)

# Wind generation constraints
wind_generation = wind_generation_constraints!(model, sets, params, wind_factors)



## Excluded constraints
#reserve_offer_curve_constraints!(model, sets, prices)
#reserve_allocation_constraints!(model, sets)
#ID_trade_limit_constraints!(model, sets, params, α = α)
#hydro_on_off_constraints!(model, sets, params)


## -- Solving model --

@info("Solving model    ("*Dates.format(now(), "HH:MM")*")")

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "TimeLimit"   => 600,
    "MIPGap" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "IntFeasTol" => 1e-9,
    "LogFile" => path*"solve_1.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)

optimize!(model)

first_solve = solve_time(model)
first_gap = relative_gap(model)
println("Solution time: ", solve_time(model), "\ngap: ", first_gap, "\n")



## -- Optimising over degenerate solutions --

@info("Optimising over degenerate solutions     ("*Dates.format(now(), "HH:MM")*")")

obj_value = objective_value(model)

# Warm start
var = all_variables(model)
var_solution = value.(var)
set_start_value.(var, var_solution)

ϵ = 0.001
obj_value_constraint = @constraint(model, objective ≥ obj_value - ϵ)

new_objective = degenerate_solutions_objective!(model, sets)
@objective(model, Min, new_objective)

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "MIPFocus"    => 3,
    "TimeLimit"   => 12000,
    "MIPGap" => 1e-6,
    "LogFile" => path*"solve_2.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)


optimize!(model)

second_solve = solve_time(model)
second_gap = relative_gap(model)
println("Solution time: ", solve_time(model), "\ngap: ", second_gap, "\n")


## -- Results summary --
@info("Saving results  ("*Dates.format(now(), "HH:MM")*")")

results_breakdown(model, sets, prices, costs, obj_value, path, run_path, α = α, first_solve = first_solve, first_gap = first_gap, second_solve = second_solve, second_gap = second_gap)

save_bid_curves(model, sets, path)

write_market_vars_jld(model, sets, path)

@info("Drawing bidding curves   ("*Dates.format(now(), "HH:MM")*")")
isdir(path*"DA_curves/") || mkdir(path*"DA_curves/")
isdir(path*"reserve_curves/") || mkdir(path*"reserve_curves/")
for t in sets.T

    plot_transposed_bid_curve(model[:x][:,t], sets.pI, t, path*"DA_curves/DA"*string(t), 
                plot_DA_prices=true, default_limits=false, price_limits =[-10, 500],
                DA_prices=prices.DA[t,:], title = "Hour "*string(t))

    if t == 23
        plot_transposed_reserve_offer_curve(model[:v][:,t], sets.pJ, t, path*"reserve_curves/reserve"*string(t),
                    plot_reserve_prices=true, default_limits=false, price_limits = [-1, 115], 
                    reserve_prices=prices.reserve[t,:], title = "Hour "*string(t))

    else
        plot_transposed_reserve_offer_curve(model[:v][:,t], sets.pJ, t, path*"reserve_curves/reserve"*string(t),
                    plot_reserve_prices=true, default_limits=false, price_limits = [-1, 25], 
                    reserve_prices=prices.reserve[t,:], title = "Hour "*string(t))
    end

end


## -- Saving optimal DA bids --

optimal_DA_bids = value.(model[:x])
optimal_DA_commitments = value.(model[:y])

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################


run_path = "DA_fixed_coordinated_"*Dates.format(now(), "yyyymmdd-HHMM")
path = path*run_path*"/"
isdir(path) || mkdir(path)
@info("path = "*path)



## -- MODEL and variables --

@info("Adding variables and objective to 2nd model")

model = Model()

market_variables!(model, sets);

hydropower_variables!(model, sets, include_start_stop = false); 

CCGT_generation_variables!(model, sets);

wind_generation_variables!(model, sets);


## -- Objective --

objective = profit_objective!(model, sets, prices, costs,
                        DA_revenue = true,
                        ID_revenue = true,
                        reserve_revenue = true,
                        imbalance_revenue =true,
                        CCGT_costs = true,
                        hydro_costs = false,
                        water_value = true)

@objective(model, Max, objective)


## -- Constraints --

@info("Building constraints")

# Market constratins
@constraint(model, fix_DA_bids, model[:x] .== optimal_DA_bids)
@constraint(model, fix_DA_commitments, model[:y] .== optimal_DA_commitments)

reserve_offer_curve = reserve_offer_curve_constraints!(model, sets, prices)

reserve_allocation = reserve_allocation_constraints!(model, sets)

capacity_limit1 = first_stage_capacity_limit_constraints!(model, sets, params)

capacity_limit2 = second_stage_capacity_limit_constraints!(model, sets, params)

α = 0.05
ID_trade_ub, ID_trade_lb = ID_trade_limit_constraints!(model, sets, params, α = α)

generation_capacity_limit = generation_capacity_constraints!(model, sets, params) 

imbalance, excess_ub, deficit_ub = imbalance_calculation_constraints!(model, sets, params)

# Hydro generation constraints

generation_conversion, discharge_lb, discharge_ub = hydro_generation_constraints!(model, sets, params) 

level_lb, level_ub, level_intertemporal = hydro_water_level_constraints!(model, sets, params)

# CCGT operation constraints
CCGT_on_off, CCGT_start, CCGT_stop = CCGT_on_off_constraints!(model, sets, params)

CCGT_min_off_time, CCGT_min_on_time = CCGT_min_on_off_time_constraints!(model, sets, params)

CCGT_generation_ub, CCGT_generation_lb = CCGT_generation_capacity_constraints!(model, sets, params)

ramping_lb, ramping_ub = CCGT_ramping_constraints!(model, sets, params)


# Wind generation constraints
wind_generation = wind_generation_constraints!(model, sets, params, wind_factors)


## Excluded constraints
# DA_offer_curve_constraints!(model, sets, prices)
# bid_quantity_order_constraints!(model, sets)
# hydro_on_off_constraints!(model, sets, params)


## -- Solving model --

@info("Solving model    ("*Dates.format(now(), "HH:MM")*")")

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "TimeLimit"   => 600,
    "MIPGap" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "IntFeasTol" => 1e-9,
    "LogFile" => path*"solve_1.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)

optimize!(model)

first_solve = solve_time(model)
first_gap = relative_gap(model)
println("Solution time: ", solve_time(model), "\ngap: ", first_gap, "\n")



## -- Optimising over degenerate solutions --

@info("Optimising over degenerate solutions     ("*Dates.format(now(), "HH:MM")*")")

obj_value = objective_value(model)

# Warm start
var = all_variables(model)
var_solution = value.(var)
set_start_value.(var, var_solution)

ϵ = 0.001
obj_value_constraint = @constraint(model, objective ≥ obj_value - ϵ)

new_objective = degenerate_solutions_objective!(model, sets)
@objective(model, Min, new_objective)

optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "MIPFocus"    => 3,
    "TimeLimit"   => 12000,
    "MIPGap" => 1e-6,
    "LogFile" => path*"solve_2.txt",
    "LogToConsole" => 0,
)
set_optimizer(model, optimizer)


optimize!(model)

second_solve = solve_time(model)
second_gap = relative_gap(model)
println("Solution time: ", solve_time(model), "\ngap: ", second_gap, "\n")


## -- Results summary --
@info("Saving results  ("*Dates.format(now(), "HH:MM")*")")

results_breakdown(model, sets, prices, costs, obj_value, path, run_path, α = α, first_solve = first_solve, first_gap = first_gap, second_solve = second_solve, second_gap = second_gap)

save_bid_curves(model, sets, path)

write_market_vars_jld(model, sets, path)
