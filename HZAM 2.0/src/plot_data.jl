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
global colors = [(:blue, 0.25), (:red, 0.25), (:purple, 0.25), (:yellow, 0.25), (:orange, 0.25), (:brown, 0.25), (:green, 0.25), (:gray, 0.25), (:cyan, 0.25), (:black, 0.25)]

global ax, points, points_active, points_inactive, mitochondria_points
global sigmoid_lines = []

function create_new_plot(hybrid_indices, mitochondria, locations)#=
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
    display(fig)=#

    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)
    
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # this sets the standard font size
    fig = Figure(resolution=(1800, 1200), figure_padding=60)
    global ax = Axis(fig[1, 1], xlabel="location", ylabel="hybrid index", title=string("HZAM simulation, generation = ", 0), xticklabelsize=45, yticklabelsize=45, titlegap=30)
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)
    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices)
    
    sigmoid_curves = []
    hybrid_zone_widths = []
    for i in 1:10
        fit = curve_fit(sigmoid, locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], initial_par)
        push!(sigmoid_curves, sigmoid(spaced_locations, fit.param))
        push!(sigmoid_lines, lines!(ax, spaced_locations, sigmoid_curves[i], color=colors[i], linewidth=20))
        push!(hybrid_zone_widths, calc_width(sigmoid_curves[i], spaced_locations))
    end

    println(hybrid_zone_widths)
    

    display(fig)
    
    #readline()
end

function update_population_plot(hybrid_indices, mitochondria, locations, generation)
    #=locations = [l.x for l in locations]
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
    ax.title = string("HZAM simulation, generatiyon = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations))=#
    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)

    delete!(ax, points)
    [delete!(ax, sigmoid_line) for sigmoid_line in sigmoid_lines]
    global sigmoid_lines = []

    global points = scatter!(ax, locations_x, locations_y, color=hybrid_indices)
    
    sigmoid_curves = []
    hybrid_zone_widths = []
    for i in 1:10
        fit = curve_fit(sigmoid, locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], initial_par)
        push!(sigmoid_curves, sigmoid(spaced_locations, fit.param))
        push!(sigmoid_lines, lines!(ax, spaced_locations, sigmoid_curves[i], color=colors[i], linewidth=20))
        push!(hybrid_zone_widths, calc_width(sigmoid_curves[i], spaced_locations))
    end

    ax.title = string("HZAM simulation, generation = ", generation)
    println("generation: ", generation, "; individuals: ", length(locations))
    println("hybrid zone width: ", sum(hybrid_zone_widths)/10)
    println("hybrid zone length: ", calc_length(sigmoid_curves, spaced_locations))
    
end


# This function produces the y value for a sigmoid, given the centre and maximum slope.
function sigmoid(x, p)  # where p[1] is cline centre, and p[2] is maxSlope 
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end


function calc_width(sigmoid_curve, locations)
    left_boundary = locations[argmin(abs.(sigmoid_curve .- 0.1))]
    right_boundary = locations[argmin(abs.(sigmoid_curve .- 0.9))]

    return right_boundary - left_boundary
end

function calc_length(sigmoid_curves, locations)
    total_length = 0.1
    mid_points = map(x -> locations[argmin(abs.(sigmoid_curves[x] .- 0.5))], collect(1:10))
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

function sort_y(y_locations)
    bins = collect(0.1:0.1:1)

    function get_indices(bin)
        return findall(y->bin-0.1 <= y < bin, y_locations)
    end

    return map(get_indices, bins)
end
end