using Plots
using CSV
using DataFrames
using Measures


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
