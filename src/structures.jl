struct Sets
    T::UnitRange{Int64}
    S::UnitRange{Int64} 
    E::UnitRange{Int64}
    W::UnitRange{Int64}
    π_S::Vector{Float64}
    π_E::Vector{Float64}
    π_W::Vector{Float64}
    I::UnitRange{Int64}
    J::UnitRange{Int64}
    pI::Vector{Float64}
    pJ::Vector{Float64}

    function Sets(;
        T = 1:24, 
        S = 1:11, 
        E = 1:6, 
        W = 1:3, 
        πS = fill(1/11, 11), 
        πE = fill(1/6, 6), 
        πW = fill(1/3, 3),
        I = 1:34,
        J = 1:17,
        pI = [-500, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 400, 4000.0],
        pJ = [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 100, 300])

        if length(πS) != S[end] || length(πE) != E[end] || length(πW) != W[end]
            throw(DomainError("Probability vector(s) have incorrect lengths."))
        end

        if length(pI) != I[end] || length(pJ) != J[end]
            throw(DomainError("Bid vector(s) have incorrect lengths."))
        end

        if pI[1] != -500 || pI[end] != 4000
            throw(DomainError("Day-ahead bid curve end points are incorrect."))
        end

        new(T, S, E, W, πS, πE, πW, I, J, pI, pJ)
    end
end

mutable struct Prices
    DA::Matrix{Float64}
    ID::Array{Float64}
    reserve::Matrix{Float64}
    up::Array{Float64}
    down::Array{Float64}

    function Prices(;
        DA = ones(2,2)*0.0,
        ID = ones(2,2)*0.0,
        reserve = ones(2,2)*0.0,
        up = ones(2,2)*0.0,
        down = ones(2,2)*0.0)

        new(DA, ID, reserve, up, down)
    end
end

struct Costs
    CCGT_startup
    CCGT_shutdown
    CCGT_generation
    Hydro_startup
    Hydro_WV                # Only single linear WV function supported
    
    function Costs(;
        CCGT_startup = 0,
        CCGT_shutdown = 0,
        CCGT_generation = 243.7,
        Hydro_startup = 0,
        Hydro_WV = [0.0015, 0])

        if !(isa(Hydro_WV, Vector{Float64}) || isa(Hydro_WV, Vector{Int64}))
            throw(DomainError("Piecewise linear water value function is currently not supported. Use single linear function"))

        elseif length(Hydro_WV) !=2
            throw(DomainError("The water value function should be defined using two parameters – its slope and intercept."))

        end

        new(CCGT_startup, CCGT_shutdown, CCGT_generation, Hydro_startup, Hydro_WV)
    end
end

struct OperationalParameters
    # Capacities
    C_CCGT              
    C_Hydro
    C_Wind
    # Hydro parameters
    eta                      # Hydro conversion coefficients
    F_max                    # Hydro minimum and maximum discharge from reservoir
    F_min
    F_inflow
    N 
    L_min                    # Hydro minimum and maximum reservoir water level
    L_max
    L_initial
    U_hydro_initial          # Hydro initial on-off status
    # CCGT parameters
    G_min_level              # CCGT generation minimum stable level
    R_CCGT                   # Maximum reserve capacity of CCGT in MW
    G_initial                # CCGT generation initial level
    A_ramp_factor            # CCGT ramping factors
    U_min_on_time            # CCGT minimum on and off times
    U_min_off_time        
    U_CCGT_initial           # CCGT initial on-off status
    dt                       # Time step length

    function OperationalParameters(;
            C_CCGT = 170,                   # Capacities        
            C_Hydro = 100.8,
            C_Wind = 40,
            eta = 0.096,                    # Hydro conversion coefficients
            F_max = 1050,                   # Hydro minimum and maximum discharge from reservoir
            F_min = 200,
            F_inflow = fill(400 , 24),      # !! Assumption nT = 24!
            N = 3600, 
            L_min = 120000000,              # Hydro minimum and maximum reservoir water level
            L_max = 136000000,
            L_initial = 128000000,
            U_hydro_initial = 1,
            G_min_level = 0.4,               # CCGT generation minimum stable level
            R_CCGT = 3,                      # Maximum reserve capacity of CCGT in MW
            G_initial = 0,                   # CCGT generation initial level
            A_ramp_factor = 1,               # CCGT ramping factors
            U_min_on_time = 7,               # CCGT minimum on and off times
            U_min_off_time = 2,        
            U_CCGT_initial = 0,              # CCGT initial on-off status
            dt = 1
    )

        new(C_CCGT, C_Hydro, C_Wind, eta, F_max, F_min, F_inflow, N, L_min, L_max, L_initial,      U_hydro_initial, G_min_level, R_CCGT, G_initial, A_ramp_factor, U_min_on_time, U_min_off_time, U_CCGT_initial, dt)
    end
end