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

    function PopulationTrackingData(
        genotypes::Vector{<:Matrix{<:Integer}},
        locations::Vector,
        male_mating_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )

        mmt_hybrid_indices = calc_traits_additive(genotypes, male_mating_trait_loci)

        sorted_indices = sort_locations_2D(locations, 50)

        @time calc_overlap(mmt_hybrid_indices, sorted_indices)

        overlap = calc_overlap(mmt_hybrid_indices, sorted_indices)

        sigmoid_curves = calc_sigmoid_curves(locations, mmt_hybrid_indices)
        mmt_cline_width = average_width(sigmoid_curves)

        population_size = length(genotypes)

        hybridity = mean(map(x -> 1 - 2 * abs(x - 0.5), mmt_hybrid_indices))

        new(population_size, hybridity, overlap, mmt_cline_width)
    end
end


"""
    OutputData

An OutputData stores all of the key data from a simulation run.

# fields
- `sim_params::SimParams`: the parameters of the simulation.
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of the population at the end of the simulation.
- `locations::Vector{<:Any}`: the locations of the population at the end of the simulation.
- `hybrid_zone_width::Real`: the cline width associated with the male mating trait
- `population_overlap::Real`: the proportion of the range occupied by males of both mating trait phenotypes
"""
struct OutputData
    sim_params::SimParams
    genotypes::Vector{<:Matrix{<:Integer}}
    locations::Vector{<:Any}
    hybrid_zone_width::Real
    population_overlap::Real
    bimodality::Real
    population_tracking_data::Vector{PopulationTrackingData}

    function OutputData(
        genotypes::Vector{<:Matrix{<:Integer}},
        locations::Vector{<:Any},
        male_mating_trait_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
        sim_params::SimParams,
        pop_track_data::Vector{PopulationTrackingData}
    )
        mmt_hybrid_indices = calc_traits_additive(genotypes, male_mating_trait_loci)

        sorted_indices = sort_locations_2D(locations, 20)

        @time calc_overlap(mmt_hybrid_indices, locations, 0.01)

        println("Overlap2: ", calc_overlap2(mmt_hybrid_indices, locations, 0.01))

        overlap = calc_overlap(mmt_hybrid_indices, locations, 0.01)

        sigmoid_curves = calc_sigmoid_curves(locations, mmt_hybrid_indices)
        mmt_cline_width = average_width(sigmoid_curves)

        locations_y = [l.y for l in locations]

        sorted_indices = sort_y(locations_y)

        bimodality = calc_bimodality_overall(
            sigmoid_curves,
            sorted_indices,
            [l.x for l in locations],
            mmt_hybrid_indices,
            0.05
        )

        new(
            sim_params,
            genotypes,
            locations,
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
For testing
"""
function calc_overlap2(
    hybrid_indices::Vector{<:Real},
    locations::Vector,
    sigma_comp::Real
)
    spaced_locations_2D = Location.(collect(0:0.01:1), collect(0:0.01:1)')

    function calc_squared_distances(focal_location, locations)
        dif_x = [l.x for l in locations] .- focal_location.x
        dif_y = [l.y for l in locations] .- focal_location.y

        return dif_x .^ 2 .+ dif_y .^ 2
    end

    function calc_density(focal_location, locations, sigma_comp)
        squared_distances = calc_squared_distances(focal_location, locations)

        return sum(exp.(-squared_distances ./ Ref(2 * (sigma_comp^2))))
    end

    species_A_locations = locations[Bool[h == 0 for h in hybrid_indices]]
    species_B_locations = locations[Bool[h == 1 for h in hybrid_indices]]

    densities_A = calc_density.(spaced_locations_2D, Ref(species_A_locations), Ref(sigma_comp))
    densities_B = calc_density.(spaced_locations_2D, Ref(species_B_locations), Ref(sigma_comp))
    total_density = calc_density.(spaced_locations_2D, Ref(locations), Ref(sigma_comp))

    min_proportion = 0.05

    return count(
        i -> densities_A[i] > min_proportion * total_density[i] &&
            densities_B[i] > min_proportion * total_density[i],
        eachindex(spaced_locations_2D)
    ) / length(spaced_locations_2D)
end

"""
function calc_overlap(
    hybrid_indices::Vector{<:Real},
    locations::Vector,
    sigma_comp::Real
)

Compute the total overlap area between the two species.

# Arguments
- `hybrid_indices::Vector{<:Real}`: the mean values of the genotypes over the loci of interest.
- `locations::Vector`: the locations of all individuals in the simulation.
- `sigma_comp::Real`: the standard deviation for the normal curve used in calculating local density.
"""
function calc_overlap(
    hybrid_indices::Vector{<:Real},
    locations::Vector,
    sigma_comp::Real
)
    spaced_locations = collect(0:0.01:1)

    function calc_overlap_in_ribbon(hybrid_indices, locations, sigma_comp)
        locations_x = [l.x for l in locations]

        function calc_density(focal_location, locations_x, sigma_comp)
            return sum(exp.(-(Ref(focal_location) .- locations_x).^2 ./ Ref(2 * (sigma_comp^2))))
        end

        species_A_locations = locations_x[Bool[h == 0 for h in hybrid_indices]]
        species_B_locations = locations_x[Bool[h == 1 for h in hybrid_indices]]

        densities_A = calc_density.(spaced_locations, Ref(species_A_locations), Ref(sigma_comp))
        densities_B = calc_density.(spaced_locations, Ref(species_B_locations), Ref(sigma_comp))
        total_density = calc_density.(spaced_locations, Ref(locations_x), Ref(sigma_comp))

        min_proportion = 0.05

        return count(
            i -> densities_A[i] > min_proportion * total_density[i] &&
                densities_B[i] > min_proportion * total_density[i],
            eachindex(spaced_locations)
        ) / length(spaced_locations)
    end

    sorted_indices = sort_y([l.y for  l in locations])

    return mean([calc_overlap_in_ribbon(hybrid_indices[i], locations[i], sigma_comp) for i in sorted_indices])
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
    sort_locations_2D(locations::Vector{<:Any}, w::Integer)

Sort the indices of the locations into a w by w matrix with each cell in the matrix 
holding a vector of indices corresponding to the individuals in that (1/w) by (1/w) portion 
of the range.
"""
function sort_locations_2D(locations::Vector{<:Any}, w::Integer)

    output = Matrix{Vector{Int64}}(undef, w, w)

    for i in eachindex(output)
        output[i] = []
    end

    for i in eachindex(locations)
        x = min(Int(floor(w * locations[i].x)) + 1, 50)
        y = min(Int(floor(w * locations[i].y)) + 1, 50)
        push!(output[x, y], i)
    end

    return output
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