"Types and functions used for analyzing data produced by the simulation."
module DataAnalysis
# export functions to be used in the PlotData module
export calc_sigmoid_curves, calc_length, calc_width, calc_bimodality_overall,
    calc_overlap_overall, spaced_locations, sort_y
export OutputData, SpatialData, PopulationTrackingData, SimParams

using LsqFit: curve_fit
using Statistics: mean

import HypothesisTests.EqualVarianceTTest

"Parameters for fitting sigmoid curves to the clines."
global initial_par = [0.0, 1.0]

"Evenly spaced locations across the range (one-dimensional)."
global spaced_locations = collect(Float32, 0:0.001:1)

"""
    SimParams

A SimParams stores the parameters used in a simulation.

# Fields
- `intrinsic_R::Real`: the intrinsic growth rate.
- `ecolDiff::Real`: the ecological difference between the two populations.
- `w_hyb::Real`: the hybrid fitness.
- `S_AM::Real`: the strength of assortative mating.
- `K_total::Integer`: the carrying capacity of the environment.
- `max_generations::Integer`: the number of generations the simulation ran for.
- `sigma_disp::Real`: the standard deviation for dispersal distance.
- `total_loci::Integer`: the number of loci in the genotypes.
- `female_mating_trait_loci`: the loci controlling mating preference.
- `male_mating_trait_loci`: the loci controlling mating cue.
- `competition_trait_loci`: the loci controlling what resources the individual uses.
- `hybrid_survival_loci`: the loci controlling hybrid fitness.
- `per_reject_cost::Real`: the loss in fitness incurred by a female after rejecting one male.
"""
struct SimParams
    intrinsic_R::Real
    ecolDiff::Real
    w_hyb::Real
    S_AM::Real
    K_total::Integer
    max_generations::Integer
    sigma_disp::Real
    total_loci::Integer
    female_mating_trait_loci
    male_mating_trait_loci
    competition_trait_loci
    hybrid_survival_loci
    per_reject_cost::Real
end

"""
    PopulationTrackingData

A PopulationTrackingData stores the size, average hybrid index, overlap, and hybrid zone 
width of the population.

As the simulation runs a vector of PopulationTrackingData is used to keep track of the key 
population data each generation.

# Fields
- `population_size::Real`: the total number of individuals.
- `hybridness::Real`: the average hybrid index in the population (0 means phenotypically pure, 1 means that half of the genetic material is from one species and half from the other).
- `overlap::Real`: the proportion of the range containing both species.
- `width::Real`: the width of the hybrid zone over all functional loci.

# Constructors
```julia
- PopulationTrackingData(
    genotypes::Vector{<:Matrix{<:Integer}},
    locations::Vector,
    loci::NamedTuple
)
```
"""
struct PopulationTrackingData
    population_size::Real
    hybridness::Real
    overlap::Base.Real
    width::Real

    function PopulationTrackingData(
        genotypes::Vector{<:Matrix{<:Integer}},
        locations::Vector,
        loci::NamedTuple
    )
        locations_x = [l.x for l in locations]
        locations_y = [l.y for l in locations]

        sorted_indices = sort_y(locations_y)

        hybrid_indices_functional = calc_traits_additive(
            genotypes,
            loci.functional
        )
        population_size = length(genotypes)
        hybridness = mean(map(x -> 1 - 2 * abs(x - 0.5), hybrid_indices_functional))
        overlap = calc_overlap_overall(locations_x, hybrid_indices_functional, sorted_indices)

        sigmoid_curves = calc_sigmoid_curves(locations, hybrid_indices_functional)

        width = average_width(sigmoid_curves)

        new(population_size, hybridness, overlap, width)
    end
end

"""
    SpatialData

A SpatialData stores the bimodality, gene flow, hybrid zone movement, cline width, and cline 
position data for the hybrid zone.

# Fields
- `bimodality::Real`: the proportion of phenotypically pure individuals within half a dispersal distance of the middle of the hybrid zone.
- `gene_flows::NamedTuple`: the proportion of genetic material originating in the other species found in individuals above a phenotypic purity cutoff.
- `cline_widths::NamedTuple`: the width of the cline for each trait.
- `cline_positions::NamedTuple`: the position of the midpoint of the cline for each trait.
- `variance::Real`: the average movement of the hybrid zone each generation.

# Constructors
```julia
- SpatialData(
        locations::Vector,
        sigma_disp::Real,
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::NamedTuple,
        last_generation::Bool
    )
- SpatialData(output_data::Vector{SpatialData})
```

# Details on behaviour of different constructors

The first constructor computes the spatial data from the population's genotypes, standard 
deviation in dispersal distance, genotypes, and which loci control which traits.

The second constructor accepts a vector of SpatialData, computes the variance in the 
functional cline position, and the average over time of all remaining fields.

"""
struct SpatialData
    bimodality::Real
    gene_flows::NamedTuple    
    cline_widths::NamedTuple    
    cline_positions::NamedTuple    
    variance::Real

    function SpatialData(
        locations::Vector,
        sigma_disp::Real,
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::NamedTuple
    )
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

        sigmoid_curves = calc_sigmoid_curves(locations, hybrid_indices_functional)

        bimodality = calc_bimodality_overall(
            sigmoid_curves,
            sorted_indices,
            locations_x,
            hybrid_indices_functional,
            sigma_disp
        )

        new(
            bimodality,
            gene_flows,
            cline_widths,
            cline_positions)
    end

    function SpatialData(output_data::Vector{SpatialData})
        bimodality = mean([o.bimodality for o in output_data])
        gene_flows = average_gene_data([o.gene_flows for o in output_data])
        positions = [o.cline_positions.functional for o in output_data]

        variance = calc_variance(positions)
        cline_widths = average_gene_data([o.cline_widths for o in output_data])
        cline_positions = average_gene_data([o.cline_positions for o in output_data])

        new(
            bimodality,
            gene_flows,
            cline_widths,
            cline_positions,
            variance
        )
    end
end

"""
    OutputData

An OutputData stores all of the key data from a simulation run.

# fields
- `sim_params::SimParams`: the parameters of the simulation.
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of the population at the end of the simulation.
- `locations::Vector{<:Any}`: the locations of the population at the end of the simulation.
- `spatial_data::Union{SpatialData,Vector{<:SpatialData}}`: the spatial characteristics of the hybrid zone over time.
- `population_tracking_data::Vector{<:PopulationTrackingData}`: the key measures of the population over time.
- `fitnesses::Vector{<:Dict}`: the average number of offspring per phenotype per generation.
- `mating_success::Vector{<:Dict}`: the average number of mates per male phenotype per generation.
- `phenotypes::Vector{NamedTuple}`: the number of individuals with each phenotype per generation.
"""
struct OutputData    
    sim_params::SimParams    
    genotypes::Vector{<:Matrix{<:Integer}}    
    locations::Vector{<:Any}    
    spatial_data::Union{SpatialData,Vector{<:SpatialData}}    
    population_tracking_data::Vector{<:PopulationTrackingData}    
    fitnesses::Vector{<:Dict}    
    mating_success::Vector{<:Dict}    
    phenotypes::Vector{NamedTuple}
end


"""
    calc_position(sigmoid_curves::Vector{<:Vector{T} where <:Real})

Calculate the x location at the middle of the cline based on a series of sigmoid curves 
fitting the cline at different ranges on the y axis.
"""
function calc_position(sigmoid_curves::Vector{<:Vector{<:Real}})
    function midpoint(sigmoid_curve)
        return spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    end

    return mean(map(midpoint, sigmoid_curves))
end

"""
    calc_variance(positions::Vector{<:Real})

Calculate the average movement of the hybrid zone each generation.
"""
function calc_variance(positions::Vector{<:Real})
    zone_movements = map(
        i -> abs(positions[i] - positions[max(1, i - 1)]),
        collect(2:length(positions))
    )
    return mean(zone_movements)
end

"""
    calc_sigmoid_curve(locations_x::Vector{<:Real}, hybrid_indices::Vector{<:Real})

Compute a sigmoid curve that models hybrid index vs location on the x axis.

Return the output of the sigmoid curve function on 1000 evenly spaced locations.
"""
function calc_sigmoid_curve(locations_x::Vector{<:Real}, hybrid_indices::Vector{<:Real})
    # calculate the parameters for the curve fitting
    fit = curve_fit(sigmoid, locations_x, hybrid_indices, initial_par)
    # compute the values of the sigmoid curve at evenly spaced x values
    return sigmoid(spaced_locations, fit.param)
end

"""
    calc_sigmoid_curves(locations::Vector, hybrid_indices::Vector{<:Real})

Divide the range into ten horizontal strips and compute a sigmoid curve fitting the hybrid 
indices along each strip.
"""
function calc_sigmoid_curves(locations::Vector, hybrid_indices::Vector{<:Real})
    locations_x = [l.x for l in locations]
    locations_y = [l.y for l in locations]

    sorted_indices = sort_y(locations_y)

    return [
        calc_sigmoid_curve(
            locations_x[sorted_indices[i]],
            hybrid_indices[sorted_indices[i]])
        for i in eachindex(sorted_indices)
    ]
end

"""
    sigmoid(x::Vector{<:Real}, p::Vector{<:Real})

Compute the y value for a sigmoid given the x value, centre, and maximum slope.

# Arguments
- `x::Real`: the value along the x axis where the sigmoid is to be calculated.
- `p::Vector{<:Real}`: the fit parameters where p[1] is the cline centre and p[2] is 
the maximum slope.
"""
function sigmoid(x::Vector{<:Real}, p::Vector{<:Real})
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

"""
    calc_width(sigmoid_curve::Vector{<:Real})

Compute the width of a cline given a sigmoid curve.

The cline width is defined as the distance between where the sigmoid curve passes through 
0.1 and 0.9.
"""
function calc_width(sigmoid_curve::Vector{<:Real})
    left_boundary = spaced_locations[argmin(abs.(sigmoid_curve .- 0.1))]
    right_boundary = spaced_locations[argmin(abs.(sigmoid_curve .- 0.9))]

    return right_boundary - left_boundary
end

"""
    average_width(sigmoid_curves::Vector{<:Vector{<:Real}})

Compute the width of a cline given a series of sigmoid curves representing the cline along 
horizontal strips of the range.
"""
function average_width(sigmoid_curves::Vector{<:Vector{<:Real}})
    mean([calc_width(sigmoid_curves[i]) for i in eachindex(sigmoid_curves)])
end

"""
    calc_length(sigmoid_curves::Vector{<:Vector{<:Real}})

Compute the approximate length of the hybrid zone.

This is done by finding the midpoints of a series of sigmoid curves representing the cline 
along evenly spaced horizontal strips of the range. The cline length is defined as the 
length of the shortest path traversing the entire range passing through each midpoint.
"""
function calc_length(sigmoid_curves::Vector{<:Vector{<:Real}})
    total_length = 0.1
    # get the x coordinates of each sigmoid curve where the curve passes through 0.5
    mid_points = map(
        x -> spaced_locations[argmin(abs.(sigmoid_curves[x] .- mean(sigmoid_curves[x])))],
        collect(1:10))
    # add up the distances between each midpoint
    for i in 2:10
        total_length += sqrt(0.1^2 + (mid_points[i] - mid_points[i-1])^2)
    end
    return total_length
end

"""
    calc_all_gene_flow(
        genotypes::Vector{<:Matrix{<:Real}},
        hybrid_indices_functional::Vector{<:Real},
        loci::NamedTuple
    )

For each trait, compute the proportion of genetic material originating in the other species 
found in individuals above a phenotypic purity cutoff.

# Arguments

- `genotypes::Vector{<:Matrix{<:Real}}`: the genotypes of all individuals.
- `hybrid_indices_functional::Vector{<:Real}`: the mean values of the genotypes over the functional loci only.
- `loci::NamedTuple`: the name and loci range of each trait of interest.
"""
function calc_all_gene_flow(
    genotypes::Vector{<:Matrix{<:Real}},
    hybrid_indices_functional::Vector{<:Real},
    loci::NamedTuple
)
    species_A_genotypes = genotypes[filter(
        x -> hybrid_indices_functional[x] < 0.25,
        eachindex(hybrid_indices_functional)
    )]

    species_B_genotypes = genotypes[filter(
        x -> hybrid_indices_functional[x] > 0.75,
        eachindex(hybrid_indices_functional)
    )]

    function calc_gene_flow(loci_range)
        species_A_indices = calc_traits_additive(species_A_genotypes, loci_range)
        species_B_indices = calc_traits_additive(species_B_genotypes, loci_range)
        (mean(species_A_indices) + mean(1 .- species_B_indices)) / 2
    end

    return map(calc_gene_flow, loci)
end

"""
    calc_all_cline_widths_and_positions(
        genotypes::Vector{<:Matrix{<:Real}},
        locations::Vector,
        loci::NamedTuple
    )

For each trait, compute the cline width and the location of the middle.

# Arguments

- `genotypes::Vector{<:Matrix{<:Real}}`: the genotypes of all individuals.
- `locations::Vector`: the locations of all individuals.
- `loci::NamedTuple`: the name and loci range of each trait of interest.
"""
function calc_all_cline_widths_and_positions(
    genotypes::Vector{<:Matrix{<:Real}},
    locations::Vector,
    loci::NamedTuple
)
    hybrid_indices_per_trait = map(t -> calc_traits_additive(genotypes, t), loci)

    sigmoid_curves_per_trait = map(
        h -> calc_sigmoid_curves(locations, h),
        hybrid_indices_per_trait
    )

    cline_widths_per_trait = map(average_width, sigmoid_curves_per_trait)

    cline_positions_per_trait = map(calc_position, sigmoid_curves_per_trait)

    return cline_widths_per_trait, cline_positions_per_trait
end

"""
    calc_overlap_in_range(
        locations_x::Vector{<:Real},
        hybrid_indices_functional::Vector{<:Real}
    )

Compute the proportion of the range containing both species.

Subdivide the range into 50 boxes and sums the area covered by boxes containing at least 
10% phenotypically pure individuals of both species.
"""
function calc_overlap_in_range(
    locations_x::Vector{<:Real},
    hybrid_indices_functional::Vector{<:Real}
)
    min_proportion = 0.1

    sorted_indices = sort_locations(locations_x, 0.02)

    """
    Compute the proportion of phenotypically pure individuals from species A.
    """
    function calc_proportion(hybrid_indices_functional, species)
        return count(x -> x == species, hybrid_indices_functional) /
               length(hybrid_indices_functional)
    end

    """
    Determine if a list of hybrid indices contains proportions of both species A and species 
    B above the cutoff.
    """
    function overlaps(hybrid_indices_functional)
        return (calc_proportion(hybrid_indices_functional, 0) > min_proportion &&
                calc_proportion(hybrid_indices_functional, 1) > min_proportion)
    end

    num_overlap_zones = count(
        x -> overlaps(hybrid_indices_functional[sorted_indices[x]]),
        eachindex(sorted_indices)
    )

    return num_overlap_zones * 0.02 * 0.1
end

"""
    function calc_overlap_overall(
        locations_x::Vector,
        hybrid_indices_functional::Vector{<:Real},
        sorted_indices::Vector{<:Vector{<:Integer}}
    )

Compute the total overlap area between the two species.

# Arguments
- `locations_x::Vector`: the x values of all the locations.
- `hybrid_indices_functional::Vector{<:Real}`: the mean values of the genotypes over the functional loci only.
- `sorted_indices::Vector{<:Vector{<:Integer}}`: the indices of the locations sorted into 
bins corresponding to non-overlapping ranges of the y values.
"""
function calc_overlap_overall(
    locations_x::Vector,
    hybrid_indices_functional::Vector{<:Real},
    sorted_indices::Vector{<:Vector{<:Integer}}
)
    sorted_locations_x = [locations_x[sorted_indices[i]] for i in eachindex(sorted_indices)]
    sorted_hybrid_indices = [
        hybrid_indices_functional[sorted_indices[i]] for i in eachindex(sorted_indices)
    ]

    overlap_per_range = calc_overlap_in_range.(sorted_locations_x, sorted_hybrid_indices)

    return sum(overlap_per_range)
end

"""
    calc_bimodality_in_range(
        sigmoid_curve::Vector{<:Real},
        locations_x::Vector{<:Real},
        hybrid_indices::Vector{<:Real},
        sigma_disp::Real
    )

Compute the proportion of phenotypically pure individuals within half a dispersal distance 
of the cline midpoint.

# Arguments
- `sigmoid_curve::Vector{<:Real}`: the sigmoid curve fitting the cline.
- `locations_x::Vector{<:Real}`: the x values of the locations.
- `hybrid_indices::Vector{<:Real}`: the mean values of the genotypes over the functional loci only.
- `sigma_disp::Real`: the standard deviation in the dispersal distance distribution.
"""
function calc_bimodality_in_range(
    sigmoid_curve::Vector{<:Real},
    locations_x::Vector{<:Real},
    hybrid_indices::Vector{<:Real},
    sigma_disp::Real
)
    center = spaced_locations[argmin(abs.(sigmoid_curve .- 0.5))]
    left = center - (sigma_disp / 2)
    right = center + (sigma_disp / 2)
    hybrid_indices_at_center =
        hybrid_indices[
            filter(i -> left <= locations_x[i] <= right, eachindex(hybrid_indices))
        ]

    bimodality = count(x -> x == 0 || x == 1, hybrid_indices_at_center) /
                 length(hybrid_indices_at_center)

    if isnan(bimodality)
        return 1
    else
        return bimodality
    end
end

"""
Compute the bimodality across the entire range (see above for further description).

    calc_bimodality_overall(
        sigmoid_curves::Vector{<:Vector{<:Real}},
        sorted_indices::Vector{<:Vector{<:Integer}},
        locations_x::Vector{<:Real},
        hybrid_indices::Vector{<:Real},
        sigma_disp::Real
    )

# Arguments
- `sigmoid_curves::Vector{<:Real}`: the sigmoid curves for each horizontal strip of the range.
- `sorted_indices::Vector{<:Vector{<:Integer}}`: the indices of the locations sorted into bins bins corresponding to non-overlapping ranges of the y values.
- `locations_x::Vector{<:Real}`: the x values of the locations.
- `hybrid_indices::Vector{<:Real}`: the mean values of the genotypes over the functional loci only.
- `sigma_disp::Real`: the standard deviation in the dispersal distance distribution.
"""
function calc_bimodality_overall(
    sigmoid_curves::Vector{<:Vector{<:Real}},
    sorted_indices::Vector{<:Vector{<:Integer}},
    locations_x::Vector{<:Real},
    hybrid_indices::Vector{<:Real},
    sigma_disp::Real
)
    sorted_locations_x = [
        locations_x[sorted_indices[i]] for i in eachindex(sorted_indices)
    ]
    sorted_hybrid_indices = [
        hybrid_indices[sorted_indices[i]] for i in eachindex(sorted_indices)
    ]

    bimodality_per_range = calc_bimodality_in_range.(
        sigmoid_curves,
        sorted_locations_x,
        sorted_hybrid_indices,
        Ref(sigma_disp)
    )

    return mean(bimodality_per_range)
end

"""
    calc_linkage_diseq(
        genotypes::Vector{<:Matrix{<:Integer}},
        l1::Integer,
        l2::Integer
    )

Compute the linkage disequilibrium (calculated using the Pearson coefficient) between two 
loci.

``r = \\frac{\\sum{(x_i -\\bar{x})(y_i -\\bar{y})}}
{\\sqrt{\\sum{(x_i -\\bar{x})^2}\\sum{(y_i -\\bar{y})^2}}}``

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotype of every individual.
- `l1::Integer`: the first focal loci.
- `l2::Integer`: the second focal loci.
"""
function calc_linkage_diseq(
    genotypes::Vector{<:Matrix{<:Integer}},
    l1::Integer,
    l2::Integer
)
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

"""
    calc_average_linkage_diseq(genotypes::Vector{<:Matrix{<:Integer}}, loci::NamedTuple)

Compute the linkage disequilibrium between each loci and return a table of values of the 
average correlation between loci of each trait.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotype of every individual.
- `loci::NamedTuple`: the name and loci range of each trait of interest.
"""
function calc_average_linkage_diseq(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple
)
    num_loci = length(loci.overall)
    n = length(loci)

    rows = (1:num_loci)
    cols = (1:num_loci)'

    linkage_diseq = calc_linkage_diseq.(Ref(genotypes), rows, cols)
    loci = [loci...]

    function average_linkage_diseq(l1, l2)
        return l1, l2, mean(linkage_diseq[loci[l1], loci[l2]])
    end

    return average_linkage_diseq.((1:n), (1:n)')
end

"""
    calc_trait_correlation(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci_range1::Union{UnitRange{<:Integer},Vector{<:Integer}},
        loci_range2::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )

Compute the correlation between two traits.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotype of every individual.
- `loci_range1::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci range of the first trait.
- `loci_range2::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci range of the second trait.
"""
function calc_trait_correlation(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci_range1::Union{UnitRange{<:Integer},Vector{<:Integer}},
    loci_range2::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    hybrid_indices1 = calc_traits_additive(genotypes, loci_range1)
    hybrid_indices2 = calc_traits_additive(genotypes, loci_range2)

    covariance = sum((hybrid_indices1 .- Ref(mean(hybrid_indices1))) .*
                     (hybrid_indices2 .- Ref(mean(hybrid_indices2))))

    deviation = sqrt(sum((hybrid_indices1 .- Ref(mean(hybrid_indices1))) .^ 2) *
                     sum((hybrid_indices2 .- Ref(mean(hybrid_indices2))) .^ 2))

    correlation = covariance / deviation

    if isnan(correlation)
        return 0
    end

    return correlation
end

"""
    calc_all_trait_correlations(genotypes::Vector{<:Matrix{<:Integer}}, loci::NamedTuple)

Compute the correlation between each trait.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotype of every individual.
- `loci::NamedTuple`: the name and loci range of each trait of interest.
"""
function calc_all_trait_correlations(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple
)
    num_traits = length(loci)
    output = []

    for (k1, v1) in collect(pairs(loci))
        for (k2, v2) in collect(pairs(loci))
            push!(output, (k1, k2, calc_trait_correlation(genotypes, v1, v2)))
        end
    end

    return output
end

"""
    sort_y(y_locations::Vector{<:Real})

Sort the y coordinate indices into 10 vectors [0, 0.1), [0.1, 0.2), etc.

# Example
```jldoctest
sort_y([0.01, 0.5, 0.24, 0.9])
10-element Vector{<:Vector{Int64}}:
 [1]
 []
 [3]
 []
 []
 [2]
 []
 []
 []
 [4]
 ```
"""
function sort_y(y_locations::Vector{<:Real})
    return sort_locations(y_locations, 0.1)
end

"""
    sort_locations(A::AbstractArray{<:Real}, bin_size::Real)

Sort the given list of numbers into bins of a given size spanning 0 to 1.
"""
function sort_locations(A::AbstractArray{<:Real}, bin_size::Real)
    bins = collect(bin_size:bin_size:1)

    function get_indices(bin)
        return findall(x -> bin - bin_size <= x < bin, A)
    end

    return map(get_indices, bins)
end

"""
    average_gene_data(gene_data::Vector{<:NamedTuple})

Compute the mean for each field of the given list of composite datatypes. 
"""
function average_gene_data(gene_data::Vector{<:NamedTuple})
    function average_trait_gene_data(trait)
        mean([gd[trait] for gd in gene_data])
    end

    traits = keys(gene_data[1])

    return (; zip(traits, map(average_trait_gene_data, traits))...)
end

"""
    calc_traits_additive(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )::Vector{Float32}

Compute the mean values of the genotypes passed to it at the given loci. Used to determine 
trait values in an additive way.
"""
function calc_traits_additive(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)::Vector{Float32}

    N = length(genotypes)
    # calculate the mean value of the genotype across the list of loci given
    function mean(genotype, loci)
        return sum(genotype[:, loci]) / (2 * length(loci))
    end

    traits = map(x -> mean(genotypes[x], loci), 1:N)
    return traits
end

"""
    average_data_per_phenotype(
        data::Vector{<:Integer},
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )

Compute the average number of mates for males of each phenotype for a given 
loci range.

# Arguments
- `data::Vector{<:Integer}`: a list of the data to be averaged for each element in the genotypes vector
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes for each individual 
- `loci::Union{UnitRange{<:Integer},Vector{<:Integer}`: the loci range governing the trait of interest.
"""
function average_data_per_phenotype(
    data::Vector{<:Integer},
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    hybrid_indices = calc_traits_additive(genotypes, loci)
    hybrid_indices = round.(hybrid_indices, digits=4)
    n = 2 * length(loci)
    different_phenotypes = Float32.(round.((0:n) .* (1 / n), digits=4))

    sum_per_phenotype = Dict{Real,Integer}(different_phenotypes .=> zeros(n + 1))

    function calc_average(i)
        return sum_per_phenotype[i] / count(x -> x == i, hybrid_indices)
    end

    for i in eachindex(data)
        sum_per_phenotype[hybrid_indices[i]] += data[i]
    end

    average_per_phenotype = calc_average.(different_phenotypes)

    return Dict(different_phenotypes .=> average_per_phenotype)
end

"""
    function count_phenotypes_at_loci(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )

Return the number of individuals with each phenotype for the given loci range.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of every individual
- `loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci range over which the phenotypes are calculated.
"""
function count_phenotypes_at_loci(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    hybrid_indices = calc_traits_additive(genotypes, loci)
    hybrid_indices = round.(hybrid_indices, digits=4)
    n = 2 * length(loci)
    different_hybrid_indices = (0:n) .* (1 / n)
    different_hybrid_indices = round.(different_hybrid_indices, digits=4)
    different_hybrid_indices = Float32.(different_hybrid_indices)

    function count_phenotype(i)
        count(x -> x == i, hybrid_indices)
    end

    return Dict(different_hybrid_indices .=> count_phenotype.(different_hybrid_indices))
end

"""
    find_fixed_alleles(genotypes::Vector{<:Matrix{<:Integer}})

Determine at which loci the population has lost its genetic diversity.
"""
function find_fixed_alleles(genotypes::Vector{<:Matrix{<:Integer}})
    extinct = []
    for i in 1:size(genotypes[1], 2)
        genotypes_at_locus = [g[:, i] for g in genotypes]
        if all(Ref(genotypes_at_locus[1]) .== genotypes_at_locus)
            push!(extinct, i)
        end
    end
    return extinct
end

"""
    calc_chi_squared(genotypes::Vector{<:Matrix{<:Integer}}, locus::Integer)
    
Compare the observed genotype frequencies at the given locus to those expected from 
Hardy-Weinberg analysis.
"""
function calc_chi_squared(genotypes::Vector{<:Matrix{<:Integer}}, locus::Integer)
    genotypes = [g[:, locus] for g in genotypes]
    alleles = vcat([g[1] for g in genotypes], [g[2] for g in genotypes])

    N = length(genotypes)
    p_A = count(x -> x == 0, alleles) / (2 * N)
    p_B = 1 - p_A
    n_AB = count(x -> x == [0; 1] || x == [1; 0], genotypes)
    n_AA = count(x -> x == [0; 0], genotypes)
    n_BB = count(x -> x == [1; 1], genotypes)


    return (((n_AA - N * p_A^2)^2) / (N * p_A^2)) +
           (((n_AB - 2 * N * p_A * p_B)^2) / (2 * N * p_A * p_B)) +
           (((n_BB - N * p_B^2)^2) / (N * p_B^2))
end
end