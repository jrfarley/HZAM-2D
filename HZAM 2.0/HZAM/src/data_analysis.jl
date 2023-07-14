module DataAnalysis
export calc_sigmoid_curve, calc_length, calc_gene_flow, calc_width, spaced_locations, calc_bimodality_overall

using LsqFit

global initial_par = [0.0, 1.0]

global spaced_locations = collect(Float32, 0:0.001:1)

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

# finds the mean value of an array
function mean(A)
    return sum(A) / length(A)
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

# calculates gene flow by finding the proportion of genes that originated in the other population in 100 individuals 
# at the far ends of the range
function calc_gene_flow(hybrid_indices, locations_x)
    species_A_indices = []
    species_B_indices = []
    for i in 1:10
        min_index = argmin(locations_x)
        deleteat!(locations_x, min_index)
        push!(species_A_indices, hybrid_indices[min_index])
        deleteat!(hybrid_indices, min_index)

        max_index = argmax(locations_x)
        deleteat!(locations_x, max_index)
        push!(species_B_indices, hybrid_indices[max_index])
        deleteat!(hybrid_indices, max_index)
    end


    gene_flow = mean(species_A_indices)
    gene_flow += mean(1 .- species_B_indices)

    return gene_flow
end

function calc_bimodality_in_range(sigmoid_curve, locations_x, hybrid_indices, sigma_disp)
    center = spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    left = center - (sigma_disp / 2)
    right = center + (sigma_disp / 2)
    hybrid_indices_at_center = hybrid_indices[filter(i -> left <= locations_x[i] <= right, eachindex(hybrid_indices))]

    count(x -> x==0 || x==1, hybrid_indices_at_center) / length(hybrid_indices_at_center)
end

function calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices, sigma_disp)
    bimodality_per_range = [calc_bimodality_in_range(sigmoid_curves[i], locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], sigma_disp) for i in eachindex(sorted_indices)]

    return mean(bimodality_per_range)
end
end