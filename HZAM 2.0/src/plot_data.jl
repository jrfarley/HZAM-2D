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

global ax, points, points_active, points_inactive, mitochondria_points, sigmoid_line

function create_new_plot(hybrid_indices_active, mitochondria_active, locations_active, hybrid_indices_inactive, mitochondria_inactive, locations_inactive)
    active_locations_x = [l.x for l in locations_active]
    active_locations_y = [l.y for l in locations_active]
    inactive_locations_x = [l.x for l in locations_inactive]
    inactive_locations_y = [l.y for l in locations_inactive]

    

    locations_x = [active_locations_x; inactive_locations_x]
    locations_y = [active_locations_y; inactive_locations_y]

    hybrid_indices = [hybrid_indices_active; hybrid_indices_inactive]

    #=fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30)
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)
    global points_active = scatter!(ax, locations_active, jitter(hybrid_indices_active), color=(:blue, 0.5))
    global points_inactive = scatter!(ax, locations_inactive, jitter(hybrid_indices_inactive), color=(:orange, 0.5))
    global mitochondria_points = scatter!(ax, [locations_active; locations_inactive], jitter([mitochondria_active; mitochondria_inactive]), color=(:green, 0.8))
    fit = curve_fit(sigmoid, [locations_active; locations_inactive], [hybrid_indices_active; hybrid_indices_inactive], initial_par)
    global sigmoid_line = lines!(ax, spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)
    display(fig)=#

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30)
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)
    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices)
    display(fig)
end

function create_new_plot(hybrid_indices, mitochondria, locations)
    locations = [l.x for l in locations]
    mitochondria = 0.25 .+ mitochondria ./ 2
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30)
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)
    global points = scatter!(ax, locations, jitter(hybrid_indices), color=(:blue, 0.5))
    global mitochondria_points = scatter!(ax, locations, jitter(mitochondria), color=(:green, 0.8))
    fit = curve_fit(sigmoid, locations, hybrid_indices, initial_par)
    global sigmoid_line = lines!(ax, spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)
    display(fig)
end

function update_population_plot(hybrid_indices_active, mitochondria_active, locations_active, hybrid_indices_inactive, mitochondria_inactive, locations_inactive, generation)#=
    # display(scatter([locations_F; locations_M], functionalLoci_HI_all_inds))
    # fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
    # lines!(spaced_locations, sigmoid(spaced_locations, fit.param))  # add the sigmoid fit to the plot
    mitochondria = 0.25 .+ [mitochondria_active; mitochondria_inactive] ./ 2
    delete!(ax, points_active)
    delete!(ax, points_inactive)
    delete!(ax, mitochondria_points)
    delete!(ax, sigmoid_line)
    global points_active = scatter!(ax, locations_active, jitter(hybrid_indices_active), color=(:blue, 0.5))
    global points_inactive = scatter!(ax, locations_inactive, jitter(hybrid_indices_inactive), color=(:orange, 0.5))
    global mitochondria_points = scatter!(ax, [locations_active; locations_inactive], jitter(mitochondria), color=(:green, 0.8))
    fit = curve_fit(sigmoid, [locations_active; locations_inactive], [hybrid_indices_active; hybrid_indices_inactive], initial_par)
    global sigmoid_line = lines!(spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)  # add the sigmoid fit to the plot
    ax.title = string("HZAM simulation, generation = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations_active)+length(locations_inactive))=#

    active_locations_x = [l.x for l in locations_active]
    active_locations_y = [l.y for l in locations_active]
    inactive_locations_x = [l.x for l in locations_inactive]
    inactive_locations_y = [l.y for l in locations_inactive]

    locations_x = [active_locations_x; inactive_locations_x]
    locations_y = [active_locations_y; inactive_locations_y]

    hybrid_indices = [hybrid_indices_active; hybrid_indices_inactive]

    delete!(ax, points)

    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices)
    ax.title = string("HZAM simulation, generation = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations_active)+length(locations_inactive))
    
end

function update_population_plot(hybrid_indices, mitochondria, locations, generation)
    locations = [l.x for l in locations]
    mitochondria = 0.25 .+ mitochondria ./ 2
    # display(scatter([locations_F; locations_M], functionalLoci_HI_all_inds))
    # fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
    # lines!(spaced_locations, sigmoid(spaced_locations, fit.param))  # add the sigmoid fit to the plot
    delete!(ax, points)
    delete!(ax, mitochondria_points)
    delete!(ax, sigmoid_line)
    global points = scatter!(ax, locations, jitter(hybrid_indices), color=(:blue, 0.5))
    global mitochondria_points = scatter!(ax, locations, jitter(mitochondria), color=(:green, 0.8))
    fit = curve_fit(sigmoid, locations, hybrid_indices, initial_par)
    global sigmoid_line = lines!(spaced_locations, sigmoid(spaced_locations, fit.param), color=(:blue, 0.25), linewidth=20)  # add the sigmoid fit to the plot
    ax.title = string("HZAM simulation, generation = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations))
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