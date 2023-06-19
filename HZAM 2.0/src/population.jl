module Population

include("plot_data.jl")

using .Plot_Data

using Test

export PopulationData, Location
export initialize_population, update_population
export choose_closest_male, calc_match_strength, generate_offspring_genotype
export plot_population, update_plot
export calc_traits_additive
export assign_zone

# stores a location with an x coordinate and a y coordinate
struct Location
    x::Float32
    y::Float32

    # creates a location the with given x and y coordinates
    function Location(x::Float32, y::Float32)
        new(x, y)
    end

    # generates a random location within the range given
    function Location(starting_range::Vector{Location})
        x = rand() * (starting_range[2].x - starting_range[1].x) + starting_range[1].x
        y = rand() * (starting_range[2].y - starting_range[1].y) + starting_range[1].y
        new(x, y)
    end

    # generates a new location
    # based on a normal distribution with width sigma_disp,
    # centred on the given location location. Constrained to be within range. 
    # geographic_limits should be a vector with two locations.
    function Location(starting_location::Location, sigma_disp::Real, geographic_limits::Vector{Location})
        x = -1
        y = -1
        while (x < geographic_limits[1].x) || (x > geographic_limits[2].x) || (y < geographic_limits[1].y) || (y > geographic_limits[2].y)
            dist = sigma_disp * randn()
            dir = rand() * 2 * pi
            x = starting_location.x + dist * cos(dir)
            y = starting_location.y + dist * sin(dir)
        end
        new(x, y)
    end
end


struct Deme
    genotypes_F::Vector{Matrix{Int8}} # the female genotypes. rows are alleles (row 1 from mother, row 2 from father) and columns are loci 
    genotypes_M::Vector{Matrix{Int8}} # the male genotypes
    locations_F::Vector{Location} # female locations
    locations_M::Vector{Location} # male locations
    mitochondria_F::Vector{Int8} # female mitochondria types (0 for species 0 and 1 for species 1)
    mitochondria_M::Vector{Int8} # male mitochondria types

    ind_useResourceA_F
    ind_useResourceA_M
    ind_useResourceB_F
    ind_useResourceB_M

    # initialize population in deme
    function Deme(starting_N,
        total_loci,
        location,
        size,
        population, ecolDiff)

        N_half = trunc(Int, starting_N)

        deme_range = [location, Location(min(location.x + size, 0.999f0), min(location.y + size, 0.999f0))]

        genotypes = fill(fill(population, 2, total_loci), N_half)

        locations = [Location(deme_range) for i in 1:N_half]

        mitochondria = fill(population, N_half)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA = calculate_ind_useResourceA.(mitochondria, ecolDiff) # at the beginning the mitochondria value is equal to the competition trait
        ind_useResourceB = calculate_ind_useResourceB.(mitochondria, ecolDiff)

        new(genotypes,
            genotypes,
            locations,
            locations,
            mitochondria,
            mitochondria,
            ind_useResourceA,
            ind_useResourceA,
            ind_useResourceB,
            ind_useResourceB)
    end

    function Deme(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M, competition_trait_loci, ecolDiff)
        # calculate new competition traits
        competition_traits_F = calc_traits_additive(genotypes_F, competition_trait_loci)
        competition_traits_M = calc_traits_additive(genotypes_M, competition_trait_loci)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA_F = calculate_ind_useResourceA.(competition_traits_F, ecolDiff)
        ind_useResourceA_M = calculate_ind_useResourceA.(competition_traits_M, ecolDiff)
        ind_useResourceB_F = calculate_ind_useResourceB.(competition_traits_F, ecolDiff)
        ind_useResourceB_M = calculate_ind_useResourceB.(competition_traits_M, ecolDiff)

        new(genotypes_F, genotypes_M, locations_F, locations_M, mitochondria_F, mitochondria_M,
            ind_useResourceA_F,
            ind_useResourceA_M,
            ind_useResourceB_F,
            ind_useResourceB_M)
    end
end


# stores all the population data (genotypes, locations, etc.) that get updated each generation
struct PopulationData
    population::Matrix{Deme}

    growth_rates_F::Matrix{Vector{Float64}} # table of the female growth rates associated with each active female index

    # initializes the genotypes, locations, mitochondria, and growth rates of the simulation
    function PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp)

        num_individuals_per_deme = K_total / 100

        intervals = collect(0.0f0:0.1f0:0.99f0)

        deme_locations = Location.(intervals, intervals')

        deme_populations = map(l -> l.x < 0.5 ? 0 : 1, deme_locations)

        demes = Deme.(Ref(num_individuals_per_deme),
            Ref(total_loci),
            deme_locations,
            Ref(0.1f0),
            deme_populations,
            Ref(ecolDiff))

        growth_rates_F = Matrix{Vector{Float64}}(undef, 10, 10)

        for zone_x in 1:10
            for zone_y in 1:10
                growth_rates_F[zone_x, zone_y] = calculate_growth_rates(demes,
                    zone_x,
                    zone_y,
                    K_total,
                    sigma_comp,
                    intrinsic_R)
            end
        end

        new(demes,
            growth_rates_F)
    end

    # initializes new collection of population data variables by replacing the old individuals with their offspring and updating the growth rates
    function PopulationData(genotypes_daughters,
        genotypes_sons,
        mitochondria_daughters,
        mitochondria_sons,
        locations_daughters,
        locations_sons,
        competition_trait_loci,
        K_total,
        sigma_comp,
        intrinsic_R,
        ecolDiff)

        demes = Deme.(genotypes_daughters,
            genotypes_sons,
            locations_daughters,
            locations_sons,
            mitochondria_daughters,
            mitochondria_sons,
            Ref(competition_trait_loci),
            Ref(ecolDiff))

        growth_rates_F = Matrix{Vector{Float64}}(undef, 10, 10)

        for zone_x in 1:10
            for zone_y in 1:10
                growth_rates_F[zone_x, zone_y] = calculate_growth_rates(demes,
                    zone_x,
                    zone_y,
                    K_total,
                    sigma_comp,
                    intrinsic_R)
            end
        end

        new(demes, growth_rates_F)
    end
end

function assign_zone(location::Location)
    zone_x = convert(Integer, trunc(10 * location.x) + 1)
    zone_y = convert(Integer, trunc(10 * location.y) + 1)
    return zone_x, zone_y
end

# calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
function calculate_ind_useResourceA(competition_traits, ecolDiff)
    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceA = competAbility .+ ((1 .- competition_traits) .* ecolDiff)

    return ind_useResourceA
end

function calculate_ind_useResourceB(competition_traits, ecolDiff)
    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceB = competAbility .+ (competition_traits .* ecolDiff)

    return ind_useResourceB
end

# calculates the ideal density (assuming normal distribution) at the location of each female
function get_ideal_densities(K_total, sigma_comp, locations_F)
    function erf(x) # approximation of the error function
        return tanh((sqrt(pi) * log(2)) * x)
    end

    function linear_density(w)
        return sqrt(pi / 2) * sigma_comp * (erf((1 - w) / (sqrt(2) * sigma_comp)) + erf(w / (sqrt(2) * sigma_comp)))
    end

    function calc_ideal_density(location) # integral from 0 to 1 of K_total*exp(-(x-focal_location)^2/(2*sigma_comp^2)) with respect to x
        return K_total * linear_density(location.x) * linear_density(location.y)
    end

    return map(calc_ideal_density, locations_F)
end

# calculates the distance between two points
function distance(location1, location2)
    return sqrt((location1.x - location2.x)^2 + (location1.y - location2.y)^2)
end

# calculates the squared distances from a set of points to a single point
function get_squared_distances(location_list, focal_location)
    dif_x = [l.x for l in location_list] .- focal_location.x
    dif_y = [l.y for l in location_list] .- focal_location.y

    return dif_x .^ 2 .+ dif_y .^ 2
end

# calculates growth rates in the spatial model
function calculate_growth_rates(population,
    zone_x,
    zone_y,
    K_total,
    sigma_comp,
    intrinsic_R)
    # in spatial model, calculate growth rates based on local resource use

    locations_F = population[zone_x, zone_y].locations_F
    neighbourhood = population[(max(zone_x - 1, 1):min(zone_x + 1, 10)), (max(zone_y - 1, 1):min(zone_y + 1, 10))]
    ind_useResourceA_all = vcat([[d.ind_useResourceB_F; d.ind_useResourceB_M] for d in neighbourhood]...)
    ind_useResourceB_all = vcat([[d.ind_useResourceB_F; d.ind_useResourceB_M] for d in neighbourhood]...)
    locations_all = vcat([[d.locations_F; d.locations_M] for d in neighbourhood]...)
    locations_F = population[zone_x, zone_y].locations_F

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    ideal_densities_at_locations_F = get_ideal_densities(K_total, sigma_comp, locations_F) # this applies the above function to each geographic location

    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

    # assume both resources have same constant density across range

    function get_useResourceA_density_real(focal_location) # this function calculates local density according to a normal curve, weighted by individual resource use
        sum(ind_useResourceA_all .* exp.(-get_squared_distances(locations_all, focal_location) ./ (2 * (sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end

    real_densities_at_locations_F_resourceA = map(get_useResourceA_density_real, locations_F) # this applies the above function to each geographic location

    function get_useResourceB_density_real(focal_location) # do the same for resource B
        sum(ind_useResourceB_all .* exp.(-get_squared_distances(locations_all, focal_location) ./ (2 * (sigma_comp^2))))
    end
    real_densities_at_locations_F_resourceB = map(get_useResourceB_density_real, locations_F) # this applies the above function to each geographic location 

    # calculate local growth rates due to each resource (according to discrete time logistic growth equation)

    local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_locations_F_resourceA ./ (ideal_densities_at_locations_F_resourceA .+ ((real_densities_at_locations_F_resourceA) .* (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_locations_F_resourceB ./ (ideal_densities_at_locations_F_resourceB .+ ((real_densities_at_locations_F_resourceB) .* (intrinsic_R - 1)))

    growth_rateA = population[zone_x, zone_y].ind_useResourceA_F .* local_growth_rates_resourceA
    growth_rateB = population[zone_x, zone_y].ind_useResourceB_F .* local_growth_rates_resourceB

    growth_rates = growth_rateA .+ growth_rateB

    return growth_rates
end

# calculates the mean value of the genotype across the list of loci given
function mean(genotype, loci)
    return sum(genotype[:, loci]) / (2 * length(loci))
end

# This function calculates each mean values of the genotypes passed to it (for each individual).
# Used to determine trait values in an additive way.
# Only those loci that are additive trait loci should be passed to this function.
function calc_traits_additive(genotypes, loci)::Vector{Float32} #=::Array{Int8,3}=#
    N = length(genotypes)
    #traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64

    traits = map(x -> mean(genotypes[x], loci), 1:N)
    return traits
end

# Finds the closest male given a list of eligible males
function choose_closest_male(elig_M, locations_M, location_mother)
    focal_male = splice!(elig_M, argmin(get_squared_distances(locations_M[elig_M], location_mother))) # this gets the index of a closest male, and removes that male from the list in elig_M
    if (distance(location_mother, locations_M[focal_male]) >= 0.1)
        return focal_male, []
    end
    return focal_male, elig_M
end

# compare male trait with female's trait (preference), and determine
# whether she accepts; note that match_strength is determined by a
# Gaussian, with a maximum of 1 and minimum of zero
function calc_match_strength(female_genotype, male_genotype, pref_SD, female_mating_trait_loci, male_mating_trait_loci)
    mating_trait_male = mean(male_genotype, male_mating_trait_loci)
    mating_trait_female = mean(female_genotype, female_mating_trait_loci)

    mating_trait_dif = mating_trait_male - mating_trait_female
    return exp((-(mating_trait_dif^2)) / (2 * (pref_SD^2)))
end

# generates the genotype for an offspring based on its parents' genotypes
function generate_offspring_genotype(mother_genotype, father_genotype)
    total_loci = size(mother_genotype, 2)
    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
    kid_genotype = Array{Int8,2}(undef, 2, total_loci)
    for locus in 1:total_loci
        kid_genotype[1, locus] = mother_genotype[rand([1 2]), locus]  # for this locus, pick a random allele from the mother
        kid_genotype[2, locus] = father_genotype[rand([1 2]), locus]  # and from the father
    end

    return kid_genotype
end


# creates a plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function plot_population(pd, functional_loci_range)
    genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
    locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]

    create_new_plot(calc_traits_additive(genotypes, functional_loci_range),
        [],
        locations)
end

# updates the plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function update_plot(pd, generation, functional_loci_range)
    genotypes = [vcat([d.genotypes_F for d in pd.population]...); vcat([d.genotypes_M for d in pd.population]...)]
    locations = [vcat([d.locations_F for d in pd.population]...); vcat([d.locations_M for d in pd.population]...)]
    update_population_plot(calc_traits_additive(genotypes, functional_loci_range),
        [],
        locations,
        generation)
end
end