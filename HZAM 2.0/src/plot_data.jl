module Plot_Data

export create_new_plot, update_population_plot

# for plotting:
#using Plots
#gr()  # use GR backend for graphs
using CategoricalArrays
using Colors, ColorSchemes
import ColorSchemes.plasma
#using Plots.PlotMeasures  # needed for plot margin adjustment
using GLMakie
using LsqFit

global geographic_limits = [0, 1]
global spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])
global functional_loci_range = 1:5
global initial_par = [0.0,1.0]

global ax, points, sigmoid_line

function create_new_plot(hybrid_indices_active, locations_active, hybrid_indices_inactive, locations_inactive)
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30)
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)
    global points_active = scatter!(ax, locations_active, jitter(hybrid_indices_active), color=(:blue, 0.5))
    global points_inactive = scatter!(ax, locations_inactive, jitter(hybrid_indices_inactive), color=(:orange, 0.5))
    fit = curve_fit(sigmoid, locations_active, hybrid_indices_active, initial_par)
    global sigmoid_line = lines!(ax, spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)
    display(fig)
    readline()
end

function update_population_plot(hybrid_indices_active, locations_active, hybrid_indices_inactive, locations_inactive, generation)
    # display(scatter([locations_F; locations_M], functionalLoci_HI_all_inds))
    # fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
    # lines!(spaced_locations, sigmoid(spaced_locations, fit.param))  # add the sigmoid fit to the plot
    delete!(ax, points_active)
    delete!(ax, points_inactive)
    delete!(ax, sigmoid_line)
    global points_active = scatter!(ax, locations_active, jitter(hybrid_indices_active), color=(:blue, 0.5))
    global points_inactive = scatter!(ax, locations_inactive, jitter(hybrid_indices_inactive), color=(:orange, 0.5))
    fit = curve_fit(sigmoid, locations_active, hybrid_indices_active, initial_par)
    global sigmoid_line = lines!(spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)  # add the sigmoid fit to the plot
    ax.title = string("HZAM simulation, generation = ", generation)
    print("generation: ")
    println(generation)
    readline()
end


# This function produces the y value for a sigmoid, given the centre and maximum slope.
function sigmoid(x, p)  # where p[1] is cline centre, and p[2] is maxSlope 
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

# This function adds jitter (small random shifts in position, to better visualize overlapping points)
function jitter(n, factor=0.02)
    n .+ (0.5 .- rand(length(n))) .* factor
end
end