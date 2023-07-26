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
using ..DataAnalysis

GLMakie.activate!(inline=false)

global colors = [(:blue, 0.25), (:red, 0.25), (:purple, 0.25), (:yellow, 0.25), (:orange, 0.25), (:brown, 0.25), (:green, 0.25), (:gray, 0.25), (:cyan, 0.25), (:black, 0.25)]

global ax, points, points_active, points_inactive, mitochondria_points
global sigmoid_lines = []
global hybrid_zone_widths = []
global hybrid_zone_positions = [0.5]

# creates the initial plot at the beginning of the simulation
function create_new_plot(hybrid_indices_all, hybrid_indices_functional, mitochondria, locations)
    locations_x = [l.x for l in locations] # x coordinates of all individuals in the simulation
    locations_y = [l.y for l in locations] # y coordinates of all individuals in the simulation

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location_x", ylabel="location_y", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels
    xlims!(-0.03, 1.03) # sets the limits for the plotted area
    ylims!(-0.03, 1.03)
    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices_functional) # adds the location of every individual to the plot

    display(fig)
end

# updates the existing plot
function update_population_plot(hybrid_indices_all, hybrid_indices_functional, mitochondria, locations, generation, sigma_disp)
    locations_x = [l.x for l in locations] # x coordinates of all individuals in the simulation
    locations_y = [l.y for l in locations] # y coordinates of all individuals in the simulation

    sorted_indices = sort_y(locations_y) # sorts the indices of the y coordinates

    delete!(ax, points) # removes the old points from the plot
    [delete!(ax, sigmoid_line) for sigmoid_line in sigmoid_lines] # removes the sigmoid curves from the plot
    global sigmoid_lines = []

    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices_functional, markersize=10) # adds the location of every individual to the plot

    
    sigmoid_curves = [calc_sigmoid_curve(locations_x[sorted_indices[i]], hybrid_indices_functional[sorted_indices[i]]) for i in eachindex(sorted_indices)]

    hybrid_zone_widths = [calc_width(sigmoid_curves[i]) for i in eachindex(sigmoid_curves)]

    [push!(sigmoid_lines, lines!(ax, spaced_locations, scale_sigmoid_curve(sigmoid_curves[i], i), color=colors[i], linewidth=20)) for i in eachindex(sigmoid_curves)]# adds the curve to the plot

    gene_flow = calc_gene_flow(hybrid_indices_all, hybrid_indices_functional)

    ax.title = string("HZAM simulation, generation = ", generation)

    
    push!(hybrid_zone_positions, calc_position(sigmoid_curves))
    
    if generation > 20
        push!(hybrid_zone_widths, sum(hybrid_zone_widths) / 10)
    end

    println("generation: ", generation, "; individuals: ", length(locations))
    println("hybrid zone width: ", sum(hybrid_zone_widths) / 10)
    println("hybrid zone length: ", calc_length(sigmoid_curves))
    println("gene flow: ", gene_flow)
    println("bimodality: ", calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices_functional, sigma_disp))
    println("overlap: ", calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices))
    println("variance: ", calc_variance(hybrid_zone_positions))
    println("")
end


# This function adds jitter (small random shifts in position, to better visualize overlapping points)
function jitter(n, factor=0.02)
    n .+ (0.5 .- rand(length(n))) .* factor
end

# scales the sigmoid curve to show up on the plot at the corresponding y-axis range
function scale_sigmoid_curve(sigmoid_curve, index)
    (sigmoid_curve .* 0.1) .+ Ref((index - 1) * 0.1)
end
end