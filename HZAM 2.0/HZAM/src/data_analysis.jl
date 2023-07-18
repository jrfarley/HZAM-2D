module DataAnalysis
export calc_sigmoid_curve, calc_length, calc_gene_flow, calc_width, calc_bimodality_overall, calc_output_data, calc_overlap_overall, calc_position, calc_variance
export average_output_data, OutputData, spaced_locations, sort_y

using LsqFit
using Statistics: mean, std

global initial_par = [0.0, 1.0]

global spaced_locations = collect(Float32, 0:0.001:1)

struct OutputData
    bimodality
    gene_flow
    hybrid_zone_width
    overlap
end

function calc_position(sigmoid_curves)
    function midpoint(sigmoid_curve)
        return spaced_locations[argmin(abs.(sigmoid_curve.-0.5))]
    end

    return mean(map(midpoint, sigmoid_curves))
end

function calc_variance(positions)
    zone_movements = map(i->abs(positions[i]-positions[max(1, i-1)]), eachindex(positions))
    return mean(zone_movements)
end

# calculates sigmoid curves that model hybrid index vs location on the x axis
# returns the output of the sigmoid curve function
function calc_sigmoid_curve(locations_x, hybrid_indices_functional)
    fit = curve_fit(sigmoid, locations_x, hybrid_indices_functional, initial_par) # calculates the parameters for the curve fitting
    return sigmoid(spaced_locations, fit.param) # gets the values of the sigmoid curve corresponding with evenly spaced x values
end


# This function produces the y value for a sigmoid, given the centre and maximum slope.
function sigmoid(x, p)  # where p[1] is cline centre, and p[2] is maxSlope 
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

# calculates the width of the hybrid zone
# the hybrid zone is defined as the region where the sigmoid curve is between 0.1 and 0.9
function calc_width(sigmoid_curve)
    left_boundary = spaced_locations[argmin(abs.(sigmoid_curve .- 0.1))]
    right_boundary = spaced_locations[argmin(abs.(sigmoid_curve .- 0.9))]

    return right_boundary - left_boundary
end

# calculates the (approximate) length of the hybrid zone
function calc_length(sigmoid_curves)
    total_length = 0.1 # 
    mid_points = map(x -> spaced_locations[argmin(abs.(sigmoid_curves[x] .- mean(sigmoid_curves[x])))], collect(1:10)) # gets the x coordinates of each sigmoid curve where the curve passes through 0.5 (the middle of the hybrid zone)
    # adds up the distance between each midpoint
    for i in 2:10
        total_length += sqrt(0.1^2 + (mid_points[i] - mid_points[i-1])^2)
    end
    return total_length
end

# calculates gene flow by finding the proportion of genes that originated in the other population among all phenotypically pure individuals
function calc_gene_flow(hybrid_indices_all, hybrid_indices_functional)
    species_A_indices = hybrid_indices_all[filter(x -> hybrid_indices_functional[x] == 0, eachindex(hybrid_indices_all))]
    species_B_indices = hybrid_indices_all[filter(x -> hybrid_indices_functional[x] == 1, eachindex(hybrid_indices_all))]

    gene_flow = (mean(species_A_indices) + mean(1 .- species_B_indices)) / 2

    return gene_flow
end

function calc_overlap_in_range(locations_x, hybrid_indices_functional)
    min_proportion = 0.1

    sorted_indices = sort_locations(locations_x, 0.02)

    function proportion_species_A(hybrid_indices_functional)
        return count(x -> x == 0, hybrid_indices_functional) / length(hybrid_indices_functional)
    end

    function proportion_species_B(hybrid_indices_functional)
        return count(x -> x == 1, hybrid_indices_functional) / length(hybrid_indices_functional)
    end

    function has_overlap(hybrid_indices_functional)
        return (proportion_species_A(hybrid_indices_functional) > min_proportion &&
                proportion_species_B(hybrid_indices_functional) > min_proportion)
    end

    num_overlap_zones = count(x->has_overlap(hybrid_indices_functional[sorted_indices[x]]), eachindex(sorted_indices))

    return num_overlap_zones * 0.02*0.1
end

function calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices)
    overlap_per_range = map(x->calc_overlap_in_range(locations_x[sorted_indices[x]], hybrid_indices_functional[sorted_indices[x]]), eachindex(sorted_indices))

    return sum(overlap_per_range)
end

function calc_bimodality_in_range(sigmoid_curve, locations_x, hybrid_indices, sigma_disp)
    center = spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    left = center - (sigma_disp / 2)
    right = center + (sigma_disp / 2)
    hybrid_indices_at_center = hybrid_indices[filter(i -> left <= locations_x[i] <= right, eachindex(hybrid_indices))]

    bimodality = count(x -> x == 0 || x == 1, hybrid_indices_at_center) / length(hybrid_indices_at_center)

    if isnan(bimodality)
        indices = filter(i -> left - 0.1 <= locations_x[i] <= right + 0.1, eachindex(locations_x))
        return 1
    else
        return bimodality
    end
end

function calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices, sigma_disp)
    bimodality_per_range = [calc_bimodality_in_range(sigmoid_curves[i], locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], sigma_disp) for i in eachindex(sorted_indices)]

    return mean(bimodality_per_range)
end

function calc_output_data(hybrid_indices_all, hybrid_indices_functional, locations, sigma_disp)
    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)

    sigmoid_curves = [calc_sigmoid_curve(locations_x[sorted_indices[i]], hybrid_indices_functional[sorted_indices[i]]) for i in eachindex(sorted_indices)]

    hybrid_zone_widths = [calc_width(sigmoid_curves[i]) for i in eachindex(sigmoid_curves)]

    gene_flows = [calc_gene_flow(hybrid_indices_all[sorted_indices[i]], hybrid_indices_functional[sorted_indices[i]]) for i in 1:10]

    bimodality = calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices_functional, sigma_disp)

    hybrid_zone_width = mean(hybrid_zone_widths)

    gene_flow = mean(gene_flows)

    overlap = calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices)

    return OutputData(bimodality, gene_flow, hybrid_zone_width, overlap)
end

# sorts the y coordinate indices into 10 vectors [0, 0.1), [0.1, 0.2), etc.
function sort_y(y_locations)
    return sort_locations(y_locations, 0.1)
end

function sort_locations(A, bin_size)
    bins = collect(bin_size:bin_size:1)

    function get_indices(bin)
        return findall(x -> bin - bin_size <= x < bin, A)
    end

    return map(get_indices, bins)
end

function average_output_data(output_data)
    bimodality = mean([o.bimodality for o in output_data])
    gene_flow = mean([o.gene_flow for o in output_data])
    hybrid_zone_width = mean([o.hybrid_zone_width for o in output_data])
    overlap = mean([o.overlap for o in output_data])

    return OutputData(bimodality, gene_flow, hybrid_zone_width, overlap)
end



end