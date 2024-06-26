## -- Cost structure --
# conventional
S_CCGT_startup = 0
S_CCGT_shutdown = 0
S_CCGT_generation = 243.7

# hydro start up cost
S_hydro_startup = 0


# Water Value function coefficients
B_hydro = [0.0015, 0]



## -- Production parameters --
# capacities
C_CCGT = 170
C_hydro = 100.8
C_wind = 40

# hydro conversion coefficients
conversion_eta = 0.096

# hydro minimum and maximum discharge from reservoir
F_max = 1050
F_min = 200
F_inflow = fill(400 , nT)
N = 3600

# hydro minimum and maximum reservoir water level
L_min = 120000000
L_max = 136000000
L_initial = 128000000

# conventional generation minimum stable level
G_min_level = 0.4

# Maximum reserve capacity of CCGT in MW
R_CCGT = 3

# conventional generation initial level
G_initial = 0

# conventional ramping factors
A_ramp_factor = 1

# conventional minimum on and off times
U_min_on_time = 7
U_min_off_time = 2

# hydro and conventional initial
U_hydro_initial = 1
U_CCGT_initial = 0


