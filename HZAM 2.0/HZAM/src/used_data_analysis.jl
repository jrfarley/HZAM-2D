"Parameters for fitting sigmoid curves to the clines."
global initial_par = [0.0, 1.0]

"Evenly spaced locations across the range (one-dimensional)."
global spaced_locations = collect(Float32, 0:0.001:1)

struct Location
    x
    y
end

"""
    SimParams

A SimParams stores the parameters used in a simulation.

# Fields
- `intrinsic_R::Real`: the intrinsic growth rate.
- `w_hyb::Real`: the hybrid fitness.
- `S_AM::Real`: the strength of assortative mating.
- `K_total::Integer`: the carrying capacity of the environment.
- `max_generations::Integer`: the number of generations the simulation ran for.
- `sigma_disp::Real`: the standard deviation for dispersal distance.
- `total_loci::Integer`: the number of loci in the genotypes.
- `female_mating_trait_loci`: the loci controlling mating preference.
- `male_mating_trait_loci`: the loci controlling mating cue.
- `hybrid_survival_loci`: the loci controlling hybrid fitness.
- `per_reject_cost::Real`: the loss in fitness incurred by a female after rejecting one male.
"""
struct SimParams
    intrinsic_R::Real
    w_hyb::Real
    S_AM::Real
    K_total::Integer
    max_generations::Integer
    sigma_disp::Real
    total_loci::Integer
    female_mating_trait_loci
    male_mating_trait_loci
    hybrid_survival_loci
    per_reject_cost::Real
end


"""
    PopulationTrackingData

A PopulationTrackingData stores the size, hybridity, overlap, and male mating trait cline 
width of the population.

As the simulation runs a vector of PopulationTrackingData is used to keep track of the key 
population data each generation.

# Fields
- `population_size::Real`: the total number of individuals.
- `hybridity::Real`: the average hybridity in the population on the male mating trait loci.
- `overlap::Real`: the proportion of the range containing both male mating trait phenotypes.
- `width::Real`: the cline width for the male mating trait.

# Constructors
```julia
- PopulationTrackingData(
    genotypes::Vector{<:Matrix{<:Integer}},
    locations::Vector,
    male_mating_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
```
"""
struct PopulationTrackingData
    population_size::Real
    hybridity::Real
    overlap::Base.Real
    width::Real
    num_offspring_A::Vector{<:Tuple}
    num_offspring_B::Vector{<:Tuple}


    function PopulationTrackingData(
        genotypes::Vector{<:Matrix{<:Integer}},
        x_locations::Vector{Float32},
        y_locations::Vector{Float32},
        male_mating_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
        overlap::Real,
        num_offspring_A::Vector{Tuple},
        num_offspring_B::Vector{Tuple}
    )
        mmt_hybrid_indices = calc_traits_additive(genotypes, male_mating_trait_loci)

        sigmoid_curves = calc_sigmoid_curves(x_locations, y_locations, mmt_hybrid_indices)
        mmt_cline_width = average_width(sigmoid_curves)

        population_size = length(genotypes)

        hybridity = mean(map(x -> 1 - 2 * abs(x - 0.5), mmt_hybrid_indices))

        new(population_size, hybridity, overlap, mmt_cline_width, num_offspring_A, num_offspring_B)
    end
end


"""
    OutputData

An OutputData stores all of the key data from a simulation run.

# fields
- `sim_params::SimParams`: the parameters of the simulation.
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of the population at the end of the simulation.
- `x_locations::Vector{Float32}`: the x coordinates of the population at the end of the simulation.
- `y_locations::Vector{Float32}`: the y coordinates of the population at the end of the simulation.
- `hybrid_zone_width::Real`: the cline width associated with the male mating trait.
- `population_overlap::Real`: the proportion of the range occupied by males of both mating trait phenotypes.
- `bimodality::Real`: the extent to which phenotypically pure individuals occur in the hybrid zone
- `population_tracking_data::Vector{PopulationTrackingData}`: the population size, hybridity, overlap, and cline width over time.
- `overlap::Real`: the proportion of the range where both species occur.
- `population_data`: the genotypes and locations for all individuals in the simulation.
"""
struct OutputData
    sim_params::SimParams
    population_data
    hybrid_zone_width::Real
    population_overlap::Real
    bimodality::Real
    population_tracking_data::Vector{PopulationTrackingData}

    function OutputData(
        genotypes::Vector{<:Matrix{<:Integer}},
        x_locations::Vector{Float32},
        y_locations::Vector{Float32},
        male_mating_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
        sim_params::SimParams,
        pop_track_data::Vector{PopulationTrackingData},
        overlap::Real,
        population_data
    )
        mmt_hybrid_indices = calc_traits_additive(genotypes, male_mating_trait_loci)

        sigmoid_curves = calc_sigmoid_curves(x_locations, y_locations, mmt_hybrid_indices)
        mmt_cline_width = average_width(sigmoid_curves)

        sorted_indices = sort_y(y_locations)

        bimodality = calc_bimodality_overall(
            sigmoid_curves,
            sorted_indices,
            x_locations,
            mmt_hybrid_indices,
            0.05
        )

        new(
            sim_params,
            population_data,
            mmt_cline_width,
            overlap,
            bimodality,
            pop_track_data
        )
    end
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
    calc_sigmoid_curves(x_locations::Vector{Float32}, y_locations::Vector{Float32}, hybrid_indices::Vector{<:Real})

Divide the range into ten horizontal strips and compute a sigmoid curve fitting the hybrid 
indices along each strip.
"""
function calc_sigmoid_curves(x_locations::Vector{Float32}, y_locations::Vector{Float32}, hybrid_indices::Vector{<:Real})
    sorted_indices = sort_y(y_locations)

    return [
        calc_sigmoid_curve(
            x_locations[sorted_indices[i]],
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
    return sort_locations(y_locations, 0.05)
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
