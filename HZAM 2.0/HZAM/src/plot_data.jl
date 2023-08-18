"Functions for plotting data while the simulation is running."
module PlotData

export create_new_plot, update_population_plot

# for plotting:
#using Plots
#gr()  # use GR backend for graphs
using Colors, ColorSchemes
import ColorSchemes.plasma
#using Plots.PlotMeasures  # needed for plot margin adjustment
using GLMakie
using ..DataAnalysis

GLMakie.activate!(inline=false) # set up the plot to display in its own window

"Colors for cline curves"
global colors = [(:blue, 0.25), (:red, 0.25), (:purple, 0.25), (:yellow, 0.25),
    (:orange, 0.25), (:brown, 0.25), (:green, 0.25), (:gray, 0.25), (:cyan, 0.25),
    (:black, 0.25)]
"The figure showing locations and hybrid indices"
global fig
"Axes on which the locations are plotted"
global ax
"Points showing the location and hybrid index of every individual"
global points
"Sigmoid curves representing the clines"
global sigmoid_lines = []

"""
    create_new_plot(
        hybrid_indices_functional::Vector,
        locations::Vector
    )

Initialize the plot of locations and hybrid indices.
"""
function create_new_plot(
    hybrid_indices_functional::Vector,
    locations::Vector,
    save_plot
)
    locations_x = [l.x for l in locations] # x coordinates of all individuals
    locations_y = [l.y for l in locations] # y coordinates of all individuals

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # set the standard font size
    global fig = Figure(resolution=(1800, 1200), figure_padding=60)
    # create the axis and labels
    global ax = Axis(
        fig[1, 1],
        xlabel="location_x",
        ylabel="location_y",
        title=string("HZAM simulation, generation = ", 0),
        xticklabelsize=45,
        yticklabelsize=45,
        titlegap=30
    )
    # set the limits for the plotted area
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)

    # add the location of every individual to the plot
    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices_functional)

    if save_plot
    dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse5")
    save(string(dir, "/", 1, ".png"), fig)
    end
    display(fig)
end

"""
    function update_population_plot(
        hybrid_indices_functional::Vector,
        locations::Vector,
        generation::Integer
    )

Update the existing plot with new locations and hybrid indices.

# Arguments
- `hybrid_indices_functional::Vector`: list of the hybrid index (value between 0 and 1) of 
every individual.
- `locations::Vector`: the location of every individual 
(must be in the same order as the hybrid indices).
- `generation::Integer`: the number of elapsed generations.
"""
function update_population_plot(
    hybrid_indices_functional::Vector,
    locations::Vector,
    generation::Integer,
    save_plot
)
    locations_x = [l.x for l in locations] # x coordinates of all individuals
    locations_y = [l.y for l in locations] # y coordinates of all individuals

    sorted_indices = sort_y(locations_y) # sort the indices of the y coordinates

    delete!(ax, points) # remove the old points from the plot
    # remove the sigmoid curves from the plot
    [delete!(ax, sigmoid_line) for sigmoid_line in sigmoid_lines]
    global sigmoid_lines = []

    # add the location of every individual to the plot
    global points = scatter!(
        ax,
        locations_x,
        locations_y,
        color=hybrid_indices_functional,
        markersize=10
    )


    sigmoid_curves = calc_sigmoid_curves(locations, hybrid_indices_functional)

    hybrid_zone_widths = [calc_width(sigmoid_curves[i]) for i in eachindex(sigmoid_curves)]

    # add the curves to the plot
    [push!(sigmoid_lines, lines!(
        ax,
        spaced_locations,
        scale_curve(sigmoid_curves[i], i),
        color=colors[i], linewidth=20
    )) for i in eachindex(sigmoid_curves)]

    ax.title = string("HZAM simulation, generation = ", generation)

    println("generation: ", generation, "; individuals: ", length(locations))
    println("hybrid zone width: ", sum(hybrid_zone_widths) / 10)
    println("hybrid zone length: ", calc_length(sigmoid_curves))
    println("bimodality: ", calc_bimodality_overall(
        sigmoid_curves,
        sorted_indices,
        locations_x,
        hybrid_indices_functional,
        0.05
    ))
    println("overlap: ", calc_overlap_overall(
        locations_x,
        hybrid_indices_functional,
        sorted_indices
    ))
    println("")

    if save_plot
    dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse5")
    save(string(dir, "/", generation, ".png"), fig)
    end
end

"""
    scale_curve(curve, index::Real)

Scale the curve output to display on the plot over the corresponding y-axis range.

# Arguments
- `curve`: the output of a function approximating the data whose range is [0,1]
- `index::Integer`: the index of the corresponding range on the y-axis. 1 refers to 
[0, 0.1), 2 to [0.1, 0.2) and so on.
"""
function scale_curve(curve, index::Integer)
    (curve .* 0.1) .+ Ref((index - 1) * 0.1)
end
end # end of Plot_Data module