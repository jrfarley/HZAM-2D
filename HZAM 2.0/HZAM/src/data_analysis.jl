module DataAnalysis
export calc_sigmoid_curves, calc_length, calc_width, calc_bimodality_overall,
    calc_overlap_overall, spaced_locations, sort_y

using LsqFit: curve_fit
using Statistics: mean, std

global initial_par = [0.0, 1.0]

global spaced_locations = collect(Float32, 0:0.001:1)

# stores all the data from a simulation
struct OutputData
    bimodality
    gene_flows
    overlap
    variance
    cline_widths
    cline_positions
    mitochondria_position
    average_linkage_diseq
end

# calculates the x location at the middle of the cline based on 
# a series of sigmoid curves of the cline at different ranges on the y axis
function calc_position(sigmoid_curves)
    function midpoint(sigmoid_curve)
        return spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    end

    return mean(map(midpoint, sigmoid_curves))
end

# calculates the average amount the hybrid zone moves each generation
function calc_variance(positions)
    zone_movements = map(i -> abs(positions[i] - positions[max(1, i - 1)]), eachindex(positions))
    return mean(zone_movements)
end

# calculates a sigmoid curve that model hybrid index vs location on the x axis
# returns the output of the sigmoid curve function on 1000 evenly spaced locations
function calc_sigmoid_curve(locations_x, hybrid_indices)
    fit = curve_fit(sigmoid, locations_x, hybrid_indices, initial_par) # calculates the parameters for the curve fitting
    return sigmoid(spaced_locations, fit.param) # gets the values of the sigmoid curve corresponding with evenly spaced x values
end

# calculates a series of sigmoid curves for 10 subsets of the range spanning 
# the x axis but evenly dividing the y axis
function calc_sigmoid_curves(locations, hybrid_indices)
    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)

    return [calc_sigmoid_curve(locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]]) for i in eachindex(sorted_indices)]
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

# calculates the width of the hybrid zone by averaging the width of the sigmoid curves that are given
function average_width(sigmoid_curves)
    mean([calc_width(sigmoid_curves[i]) for i in eachindex(sigmoid_curves)])
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
function calc_all_gene_flow(genotypes, hybrid_indices_functional, loci)
    species_A_genotypes = genotypes[filter(x -> hybrid_indices_functional[x] < 0.25, eachindex(hybrid_indices_functional))]
    species_B_genotypes = genotypes[filter(x -> hybrid_indices_functional[x] > 0.75, eachindex(hybrid_indices_functional))]

    function calc_gene_flow(loci_range)
        species_A_indices = calc_traits_additive(species_A_genotypes, loci_range)
        species_B_indices = calc_traits_additive(species_B_genotypes, loci_range)
        (mean(species_A_indices) + mean(1 .- species_B_indices)) / 2
    end

    return map(calc_gene_flow, loci)
end

# calculates gene flow by finding the proportion of genes that originated in the other population among all phenotypically pure individuals
function calc_all_cline_widths_and_positions(genotypes, locations, loci)
    hybrid_indices_per_trait = map(t -> calc_traits_additive(genotypes, t), loci)

    sigmoid_curves_per_trait = map(
        h -> calc_sigmoid_curves(locations, h),
        hybrid_indices_per_trait
    )

    cline_widths_per_trait = map(average_width, sigmoid_curves_per_trait)

    cline_positions_per_trait = map(calc_position, sigmoid_curves_per_trait)

    return cline_widths_per_trait, cline_positions_per_trait
end

# calculates the overlap in a range by subdividing the range into 50 boxes along the x axis and 
# finding the proportion of boxes that contain at least 10% "pure" individuals from both species
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

    num_overlap_zones = count(x -> has_overlap(hybrid_indices_functional[sorted_indices[x]]), eachindex(sorted_indices))

    return num_overlap_zones * 0.02 * 0.1
end

# calculates the total overlap area between the two species
function calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices)
    overlap_per_range = map(x -> calc_overlap_in_range(locations_x[sorted_indices[x]], hybrid_indices_functional[sorted_indices[x]]), eachindex(sorted_indices))

    return sum(overlap_per_range)
end

# calculates the bimodality in the range
# bimodality is defined here by the proportion of phenotypically pure individuals within half a dispersal distance of the middle of the hybrid zone
function calc_bimodality_in_range(sigmoid_curve, locations_x, hybrid_indices, sigma_disp)
    center = spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    left = center - (sigma_disp / 2)
    right = center + (sigma_disp / 2)
    hybrid_indices_at_center = hybrid_indices[filter(i -> left <= locations_x[i] <= right, eachindex(hybrid_indices))]

    bimodality = count(x -> x == 0 || x == 1, hybrid_indices_at_center) / length(hybrid_indices_at_center)

    if isnan(bimodality)
        return 1
    else
        return bimodality
    end
end

# calculates the total bimodality in the range
function calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices, sigma_disp)
    bimodality_per_range = [calc_bimodality_in_range(sigmoid_curves[i], locations_x[sorted_indices[i]], hybrid_indices[sorted_indices[i]], sigma_disp) for i in eachindex(sorted_indices)]

    return mean(bimodality_per_range)
end

# calculates the Pearson coefficient between two loci
function calc_linkage_diseq(genotypes, l1, l2)
    if l2 != l1
        genotypes = [g[:, [l1, l2]] for g in genotypes]

        haplotypes = vcat([g[1, :] for g in genotypes], [g[2, :] for g in genotypes])

        p_A = count(h -> h[1] == 0, haplotypes) / length(haplotypes)
        p_B = count(h -> h[2] == 0, haplotypes) / length(haplotypes)

        p_AB = count(h -> h == [0, 0], haplotypes) / length(haplotypes)
        D = (p_AB - (p_A * p_B))
        pearson_coefficient = (D^2) / (p_A * (1 - p_A) * p_B * (1 - p_B))
        return pearson_coefficient
    else
        return 1
    end
end

# calculates the linkage disequilibrium between loci using Pearson coefficients
# and returns a table of values representing the average correlation between two traits
function calc_average_linkage_diseq(genotypes, loci)
    num_loci = length(loci.overall)

    rows = (1:num_loci)
    cols = (1:num_loci)'

    linkage_diseq = calc_linkage_diseq.(Ref(genotypes), rows, cols)
    loci = [loci...]

    function average_linkage_diseq(l1, l2)
        return l1, l2, mean(linkage_diseq[loci[l1], loci[l2]])
    end

    return average_linkage_diseq.((1:7), (1:7)')
end

# calculates the key statistics to determine what's happening in the simulation
function calc_output_data(locations, sigma_disp, genotypes, mitochondria, loci, last_generation)
    average_linkage_diseq = missing
    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)

    hybrid_indices_functional = calc_traits_additive(
        genotypes,
        loci.functional
    )

    gene_flows = calc_all_gene_flow(
        genotypes,
        hybrid_indices_functional,
        loci
    )

    cline_widths, cline_positions = calc_all_cline_widths_and_positions(
        genotypes,
        locations,
        loci
    )

    if last_generation
        average_linkage_diseq = calc_average_linkage_diseq(
            genotypes,
            loci
        )
    end

    sigmoid_curves = calc_sigmoid_curves(locations, hybrid_indices_functional)

    bimodality = calc_bimodality_overall(sigmoid_curves, sorted_indices, locations_x, hybrid_indices_functional, sigma_disp)

    overlap = calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices)

    mitochondria_position = calc_position(calc_sigmoid_curves(locations, mitochondria))

    return OutputData(
        bimodality,
        gene_flows,
        overlap,
        missing,
        cline_widths,
        cline_positions,
        mitochondria_position,
        average_linkage_diseq)
end

# sorts the y coordinate indices into 10 vectors [0, 0.1), [0.1, 0.2), etc.
function sort_y(y_locations)
    return sort_locations(y_locations, 0.1)
end

# sorts a given list of numbers into bins of a given size spanning 0 to 1
function sort_locations(A, bin_size)
    bins = collect(bin_size:bin_size:1)

    function get_indices(bin)
        return findall(x -> bin - bin_size <= x < bin, A)
    end

    return map(get_indices, bins)
end

# calculates the average gene flow over the last 20 generations of the simulation for each trait
function average_gene_data(gene_data)
    function average_trait_gene_data(trait)
        mean([gd[trait] for gd in gene_data])
    end

    return map(average_trait_gene_data, keys(gene_data[1]))
end

# averages the output data over the last 20 generations of the simulation
function average_output_data(output_data)
    last_output = last(output_data)
    bimodality = mean([o.bimodality for o in output_data])
    gene_flows = average_gene_data([o.gene_flows for o in output_data])
    overlap = mean([o.overlap for o in output_data])
    positions = [o.cline_positions.functional for o in output_data]

    variance = calc_variance(positions)
    cline_widths = average_gene_data([o.cline_widths for o in output_data])
    cline_positions = average_gene_data([o.cline_positions for o in output_data])
    mitochondria_position = mean([o.mitochondria_position for o in output_data])

    average_linkage_diseq = last_output.average_linkage_diseq

    return OutputData(
        bimodality,
        gene_flows,
        overlap,
        variance,
        cline_widths,
        cline_positions,
        mitochondria_position,
        average_linkage_diseq
    )
end

"""
    calc_traits_additive(genotypes::Vector, loci)::Vector{Float32}

Compute the mean values of the genotypes passed to it at the given loci. Used to determine 
trait values in an additive way.
"""
function calc_traits_additive(genotypes, loci)::Vector{Float32} #=::Array{Int8,3}=#
    N = length(genotypes)
    #traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64

    # calculates the mean value of the genotype across the list of loci given
    function mean(genotype, loci)
        return sum(genotype[:, loci]) / (2 * length(loci))
    end

    traits = map(x -> mean(genotypes[x], loci), 1:N)
    return traits
end
end