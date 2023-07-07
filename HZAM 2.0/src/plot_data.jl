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

global spaced_locations = collect(Float32, 0:0.001:1)
global initial_par = [0.0,1.0]
global colors = [(:blue, 0.25), (:red, 0.25), (:purple, 0.25), (:yellow, 0.25), (:orange, 0.25), (:brown, 0.25), (:green, 0.25), (:gray, 0.25), (:cyan, 0.25), (:black, 0.25)]

global ax, points, points_active, points_inactive, mitochondria_points
global sigmoid_lines = []
global hybrid_zone_widths = []

# creates the initial plot at the beginning of the simulation
function create_new_plot(hybrid_indices, mitochondria, locations)
    locations_x = [l.x for l in locations] # x coordinates of all individuals in the simulation
    locations_y = [l.y for l in locations] # y coordinates of all individuals in the simulation

    sorted_indices = sort_y(locations_y) # sorts the indices of the y coordinates so that sorted_indices[1] is all the indices matched with y coordinates between 0 and 0.1
    # sorted indices[2] corresponds to y coordinates between 0.1 and 0.2 and so on 
    
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30) # creates the axes and labels
    xlims!(-0.03, 1.03) # sets the limits for the plotted area
    ylims!(-0.03, 1.03)
    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices) # adds the location of every individual to the plot
    
    sigmoid_curves = [] 
    hybrid_zone_widths = [] # stores the width of the hybrid zone for each interval on the y Axis
    # width is determined by the distance between where the sigmoid curve passes through 0.1 and where it passes through 0.9
    for i in 1:10
        fit = curve_fit(sigmoid, locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], initial_par) # calculates the parameters for the curve fitting
        push!(sigmoid_curves, sigmoid(spaced_locations, fit.param)) # gets the values of the sigmoid curve corresponding with evenly spaced x values
        push!(sigmoid_lines, lines!(ax, spaced_locations, sigmoid_curves[i], color=colors[i], linewidth=20)) # adds the curve to the plot
        push!(hybrid_zone_widths, calc_width(sigmoid_curves[i], spaced_locations)) # calculates the width of the hybrid zone
    end
    

    display(fig)
end

# updates the existing plot
function update_population_plot(hybrid_indices, mitochondria, locations, generation)
    locations_x = [l.x for l in locations] # x coordinates of all individuals in the simulation
    locations_y = [l.y for l in locations] # y coordinates of all individuals in the simulation

    sorted_indices = sort_y(locations_y) # sorts the indices of the y coordinates

    delete!(ax, points) # removes the old points from the plot
    [delete!(ax, sigmoid_line) for sigmoid_line in sigmoid_lines] # removes the sigmoid curves from the plot
    global sigmoid_lines = []

    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices) # adds the location of every individual to the plot
    
    sigmoid_curves = []
    hybrid_zone_widths = [] # stores the width of the hybrid zone for each interval on the y Axis
    for i in 1:10
        fit = curve_fit(sigmoid, locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], initial_par) # calculates the parameters for the curve fitting
        push!(sigmoid_curves, sigmoid(spaced_locations, fit.param)) # gets the values of the sigmoid curve corresponding with evenly spaced x values
        push!(sigmoid_lines, lines!(ax, spaced_locations, sigmoid_curves[i], color=colors[i], linewidth=20)) # adds the curve to the plot
        push!(hybrid_zone_widths, calc_width(sigmoid_curves[i], spaced_locations)) # calculates the width of the hybrid zone
    end

    ax.title = string("HZAM simulation, generation = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations))
    println("hybrid zone width: ", sum(hybrid_zone_widths)/10)
    println("hybrid zone length: ", calc_length(sigmoid_curves, spaced_locations))
    if generation > 20
        push!(hybrid_zone_widths, sum(hybrid_zone_widths)/10)
    end
    if generation==200
        println("#########################################")
        println(sum(hybrid_zone_widths)/length(hybrid_zone_widths))
        println("########################################")
    end
end   


# This function produces the y value for a sigmoid, given the centre and maximum slope.
function sigmoid(x, p)  # where p[1] is cline centre, and p[2] is maxSlope 
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

# calculates the width of the hybrid zone
# the hybrid zone is defined as the region where the sigmoid curve is between 0.1 and 0.9
function calc_width(sigmoid_curve, locations)
    left_boundary = locations[argmin(abs.(sigmoid_curve .- 0.1))]
    right_boundary = locations[argmin(abs.(sigmoid_curve .- 0.9))]

    return right_boundary - left_boundary
end

# calculates the (approximate) length of the hybrid zone
function calc_length(sigmoid_curves, locations)
    total_length = 0.1 # 
    mid_points = map(x -> locations[argmin(abs.(sigmoid_curves[x] .- 0.5))], collect(1:10)) # gets the x coordinates of each sigmoid curve where the curve passes through 0.5 (the middle of the hybrid zone)
    # adds up the distance between each midpoint
    for i in 2:10
        total_length += sqrt(0.1^2 + (mid_points[i]-mid_points[i-1])^2)
    end
    return total_length
end

function get_sigmoid_curve(locations_x, hybrid_indices)
    fit = curve_fit(sigmoid, locations_x, hybrid_indices, initial_par)
    return sigmoid(spaced_locations, fit.param)
end

# This function adds jitter (small random shifts in position, to better visualize overlapping points)
function jitter(n, factor=0.02)
    n .+ (0.5 .- rand(length(n))) .* factor
end

# sorts the y coordinate indices into 10 vectors [0, 0.1), [0.1, 0.2), etc.
function sort_y(y_locations)
    bins = collect(0.1:0.1:1)

    function get_indices(bin)
        return findall(y->bin-0.1 <= y < bin, y_locations)
    end

    return map(get_indices, bins)
end
end