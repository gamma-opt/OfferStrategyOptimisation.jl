using Plots
using CSV
using DataFrames
using Measures


function plot_bid_curve(x, pI, path::String; plot_DA_prices::Bool=false, default_limits::Bool=true, DA_prices=[], ylimits = [])

    bids = [value.(x)...]

    if default_limits
        plot(bids, pI, label = "", c = :orange, lw = 1.5, ylim = (-10, 200), xlim = (0, 350), legend=:bottomright,  dpi=200, size = (600, 500))
        xticks!(25*[1:14]...)
    else
        plot(bids, pI, label = "", c = :orange, lw = 1.5, ylim = ylimits, legend=:bottomright,  dpi=200, size = (600, 500))
    end

    if plot_DA_prices
        hline!(DA_prices, c = :grey, lw = 0.5, label = "") 
    end

    scatter!(bids, pI, c = :red, label="")
    xlabel!("MWh")
    ylabel!("Euros / MWh")
    
    savefig(path)

end



function plot_transposed_bid_curve(x, pI, h::Int, path::String; plot_DA_prices::Bool=false, default_limits::Bool=true, DA_prices=[], price_limits = [], title::String="")

    bids = [value.(x)...]

    if default_limits
        plot(pI, bids, label = "", c = :orange, lw = 1.5, ylim = (0, 350), xlim = (-10, 200), legend=:topleft,  dpi=200, size = (600, 500))
        xticks!(25*[1:14]...)
    else
        plot(pI, bids, label = "", c = :orange, lw = 1.5, xlim = price_limits, legend=:topleft,  dpi=200, size = (600, 500))
    end

    if plot_DA_prices
        vline!(DA_prices, c = :grey, lw = 0.5, label = "DA price scenarios") 
    end

    hline!([310.8], c=:black, lw=1, label="Max capacity")
    scatter!(pI, bids, c = :red, label="")
    ylabel!("MWh")
    xlabel!("Euros / MWh")
    if title == ""
        title!("Hour "*string(h))
    else
        title!(title)
    end

    
    savefig(path)

end


function plot_hydro_water_levels(l_hydro, T, S, E, path::String; title::String="Water level in reservoir", ylimits = [1.15e8, 1.4e8])

    time = [T...]

    plot(ylim = ylimits, xlim = (0, 25), legend=:bottomleft,  dpi=200, size = (600, 500))
    hline!([1.2e8, 1.36e8], label = "Water reservoir limits", c=:black, lw= 0.5)
    for s in S, e in E
        levels = [value.(l_hydro[:,s,e])...]

        plot!(time, levels, lw = 1.5, label = "")

    end
    
    xlabel!("time")
    ylabel!("Stored water (m^3)")
    title!(title)
    savefig(path)
end

function plot_transposed_reserve_offer_curve(v, pJ, h::Int, path::String; plot_reserve_prices::Bool=false, default_limits::Bool=true, reserve_prices=[], price_limits = [0, 200], title::String="")

    bids = [value.(v)...]
    curve = cumsum(bids)

    if default_limits
        plot(pJ, curve, label = "", linetype=:steppost, c = :blue, lw = 1.5, ylim = (0, 104), xlim = (-10, 200), legend=:topleft,  dpi=200, size = (600, 500))
        xticks!(25*[1:14]...)
    else
        plot(pJ, curve, label = "", linetype=:steppost, c = :blue, lw = 1.5, xlim = price_limits, legend=:topleft,  dpi=200, size = (600, 500))
    end

    if plot_reserve_prices
        vline!(reserve_prices, c = :grey, lw = 0.5, label = "Reserve price scenarios") 
    end

    hline!([100.8+3], c=:black, lw=1, label="Max reserve capacity")
    scatter!(pJ, curve, c = :blue, label="")
    ylabel!("MW")
    xlabel!("Euros / MW")
    if title == ""
        title!("Hour "*string(h))
    else
        title!(title)
    end

    
    savefig(path)

end


function plot_generation_mix(path::String, run_path::String, identifier::String; ylimit::Int = 170, legend::Bool=true, save_on_path = true)

    generation = CSV.read(path*"Generation_"*run_path*".csv", DataFrame)

    hydro = generation.Hydro

    hydro_wind = generation.Hydro .+ generation.Wind

    hydro_wind_CCGT = generation.Hydro .+ generation.Wind .+ generation.CCGT

    DA_commitment = generation.DA_commitment

    plot([0, generation.Hour...], [0, hydro_wind_CCGT...],  linetype=:steppre, fillrange = [0, 0], c = :darkgoldenrod3, label = "CCGT")
    plot!([0, generation.Hour...], [0, hydro_wind...],      linetype=:steppre, fillrange = [0, 0], c = :darkseagreen3, label = "Wind",)
    plot!([0, generation.Hour...], [0, hydro...],           linetype=:steppre, fillrange = [0, 0], c = :mediumblue, label = "Hydro",)
    plot!([0, generation.Hour...], [0, DA_commitment...],   linetype=:steppre, c = :black, lw = 2, label = "DA commitment")

    plot!(xlim = (0, 24), ylim = (0, ylimit), dpi=200, margin = 0.3cm, guidefontsize =  16, tickfontsize= 12, tick_direction = :none)

    xticks!([1:24...] .- 0.5, [string(a) for a in 1:24])
    yticks!([0:10...].*20)
    xlabel!("Hour")
    ylabel!("Generation (MWh)")

    if !legend
        plot!(legend_position = :none, size = (800, 400))
    else
        plot!(legend_position = :outertopright, legend_font_pointsize = 16, size = (1200, 450))
    end

    if save_on_path
        savefig(path*"Generation_plot_"*identifier)
    else
        savefig("Generation_plot_"*identifier)
    end
end


function plot_reserve_mix(path::String, run_path::String, identifier::String; ylimit::Int = 80, legend::Bool=true, save_on_path = true)

    generation = CSV.read(path*"Generation_"*run_path*".csv", DataFrame)

    hydro = generation.Reserve_hydro

    hydro_CCGT = generation.Reserve_hydro .+ generation.Reserve_CCGT

    Reserve_commitment = generation.Reserve_commitment

    plot([0, generation.Hour...], [0, hydro_CCGT...], linetype=:steppre, fillrange = [0, 0], c = :darkgoldenrod3, label = "CCGT")
    plot!([0, generation.Hour...], [0, hydro...], linetype=:steppre, fillrange = [0, 0], c = :steelblue, label = "Hydro")
    plot!([0, generation.Hour...], [0, Reserve_commitment...], linetype=:steppre, c = :grey18, lw = 2, label = "Reserve commitment")
    
    plot!(xlim = (0, 24), ylim = (0, ylimit), dpi=200, margin = 0.3cm, guidefontsize =  16, tickfontsize= 12, tick_direction = :none)
    
    xticks!([1:24...] .- 0.5, [string(a) for a in 1:24])
    yticks!([0:10...].*10)
    xlabel!("Hour")
    ylabel!("Reserve (MWh)")

    if !legend
        plot!(legend_position = :none, size = (800, 400))
    else
        plot!(legend_position = :outertopright, legend_font_pointsize = 16, size = (1200, 450))
    end

    if save_on_path
        savefig(path*"Reserve_plot_"*identifier)
    else
        savefig("Reserve_plot_"*identifier)
    end
end