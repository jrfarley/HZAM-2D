module Population

include("plot_data.jl")

using .Plot_Data

using Test

export PopulationData, Location
export initialize_population, update_population
export choose_closest_male, calc_match_strength, generate_offspring_genotype, disperse_individual
export plot_population, update_plot
export calc_traits_additive

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
        while (x < geographic_limits[1].x) || (x > geographic_limits[2].x)
            x = starting_location.x + (sigma_disp * randn())
        end
        while (y < geographic_limits[1].y) || (y > geographic_limits[2].y)
            y = starting_location.y + (sigma_disp * randn())
        end
        new(x, y)
    end
end


# stores all the population data (genotypes, locations, etc.) that get updated each generation
struct PopulationData
    genotypes_F::Vector{Matrix{Int8}} # the female genotypes. rows are alleles (row 1 from mother, row 2 from father) and columns are loci 
    genotypes_M::Vector{Matrix{Int8}} # the male genotypes
    locations_F::Vector{Location} # female locations
    locations_M::Vector{Location} # male locations
    mitochondria_F::Vector{Int8} # female mitochondria types (0 for species 0 and 1 for species 1)
    mitochondria_M::Vector{Int8} # male mitochondria types

    growth_rates_F::Dict{Int64,Float64} # table of the female growth rates associated with each active female index

    individuals_per_zone_F::Matrix{Vector{Int32}} # indices of the females in each zone
    individuals_per_zone_M::Matrix{Vector{Int32}} # indices of the males in each zone

    left_boundary::Int8 # zone number of the furthest left active zone
    right_boundary::Int8 # zone number of the furthest right active zone

    bottom_boundary::Int8
    top_boundary::Int8

    active_F::Vector{Int32} # indices of the currently active females (in the genotypes_F, locations_F, and mitochondria_F arrays)
    active_M::Vector{Int32} # indices of the currently active males
    inactive_F::Vector{Int32} # indices of the currently inactive females (in the inactive_genotypes_F, inactive_locations_F, and inactive_mitochondria_F arrays)
    inactive_M::Vector{Int32} # indices of the currently inactive males
    buffer_M::Vector{Vector{Int32}} # indices of the males in the buffer zone between active and inactive regions (in the inactive arrays)

    # data that only gets updated every 10 generations (used to store information about inactive individuals)
    inactive_locations_F::Vector{Location} # inactive female locations
    inactive_locations_M::Vector{Location} # inactive male locations
    inactive_genotypes_F::Vector{Matrix{Int8}} # inactive female genotypes
    inactive_genotypes_M::Vector{Matrix{Int8}} # inactive male genotypes
    inactive_mitochondria_F::Vector{Int8} # inactive female mitochondria types
    inactive_mitochondria_M::Vector{Int8} # inactive male mitochondria types

    # initializes the genotypes, locations, mitochondria, and growth rates of the simulation
    function PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp,
        geographic_limits,
        starting_range_pop0,
        starting_range_pop1,
        optimize)

        # specify ecological resource competitive abilities for two resources A and B 
        # ecolDiff = 1.0 # this is "E" in the paper 
        competAbility_useResourceA_pop0 = (1 + ecolDiff) / 2    # equals 1 when ecolDiff = 1   
        competAbility_useResourceB_pop0 = 1 - competAbility_useResourceA_pop0
        competAbility_useResourceA_pop1 = competAbility_useResourceB_pop0   # equals 0 when ecolDiff = 1
        competAbility_useResourceB_pop1 = competAbility_useResourceA_pop0

        K_A = convert(Int64, K_total / 2)  # EVEN NUMBER; carrying capacity (on resource alpha) of entire range (for two sexes combined), regardless of species 
        K_B = K_A   # EVEN NUMBER; carrying capacity (on resource beta) of entire range (for two sexes combined), regardless of species


        pop0_starting_area = (starting_range_pop0[2].x - starting_range_pop0[1].x) * (starting_range_pop0[2].y - starting_range_pop0[1].y)
        pop1_starting_area = (starting_range_pop1[2].x - starting_range_pop1[1].x) * (starting_range_pop1[2].y - starting_range_pop1[1].y)
        total_area = (geographic_limits[2].x - geographic_limits[1].x) * (geographic_limits[2].y - geographic_limits[1].y)
        pop0_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop0) + (K_B * competAbility_useResourceB_pop0)) * pop0_starting_area / total_area)
        pop0_starting_N_half = Int(round(pop0_starting_N / 2))  # starting number for each sex
        pop1_starting_N = round((2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop1) + (K_B * competAbility_useResourceB_pop1)) * pop0_starting_area / total_area)
        pop1_starting_N_half = Int(round(pop1_starting_N / 2))

        # Generate genotype array for population of females:
        # this is a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
        # and columns (D2) are loci
        # functional loci are first, followed by neutral location_index
        genotypes_F = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
        genotypes_M = genotypes_F

        # generate mitochondria array (0 for species 0 and 1 for species 1)
        mitochondria_F = [fill(0, pop0_starting_N_half); fill(1, pop1_starting_N_half)]
        mitochondria_M = mitochondria_F

        # Generate breeding locations of individuals
        locations_F, locations_M = generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA_F = [fill(competAbility_useResourceA_pop0, pop0_starting_N_half); fill(competAbility_useResourceA_pop1, pop1_starting_N_half)]
        ind_useResourceA_M = ind_useResourceA_F
        ind_useResourceB_F = [fill(competAbility_useResourceB_pop0, pop0_starting_N_half); fill(competAbility_useResourceB_pop1, pop1_starting_N_half)]
        ind_useResourceB_M = ind_useResourceB_F



        if !optimize # calculate growth rates of all individuals and initialize population without the active/inactive-zone-related variables
            growth_rates_F = calculate_growth_rates(ind_useResourceA_F,
                ind_useResourceB_F,
                ind_useResourceA_M,
                ind_useResourceB_M,
                locations_F,
                locations_M,
                K_total,
                sigma_comp,
                intrinsic_R,
                1:length(locations_F))

            new(genotypes_F,
                genotypes_M,
                locations_F,
                locations_M,
                mitochondria_F,
                mitochondria_M,
                growth_rates_F)
        else
            # group individuals into zones
            # individuals in the middle zones (5 and 6) become active whereas other individuals are inactive
            # males in zones adjacent to active zones become part of the buffer which is treated as active for mating, but remains fixed from generation to generation
            individuals_per_zone_F, individuals_per_zone_M = assign_zones(locations_F, locations_M)

            active_F = vcat(individuals_per_zone_F[5:6, :]...)
            inactive_F = vcat(individuals_per_zone_F[[(1:4); (7:10)], :]...)
            buffer_M = [vcat(individuals_per_zone_M[4, :]...), vcat(individuals_per_zone_M[7, :]...), [], []]
            active_M = vcat(individuals_per_zone_M[5:6, :]...)
            active_M = vcat(active_M, buffer_M[1], buffer_M[2])
            inactive_M = vcat(individuals_per_zone_M[[(1:3); (8:10)], :]...)
            inactive_locations_F, inactive_locations_M = (copy(locations_F), copy(locations_M))

            inactive_genotypes_F, inactive_genotypes_M = (copy(genotypes_F), copy(genotypes_M))
            inactive_mitochondria_F, inactive_mitochondria_M = (copy(mitochondria_F), copy(mitochondria_M))

            # calculate growth rates of active females only
            growth_rates_F = calculate_growth_rates(ind_useResourceA_F,
                ind_useResourceB_F,
                ind_useResourceA_M,
                ind_useResourceB_M,
                locations_F,
                locations_M,
                K_total,
                sigma_comp,
                intrinsic_R,
                active_F)

            new(genotypes_F,
                genotypes_M,
                locations_F,
                locations_M,
                mitochondria_F,
                mitochondria_M,
                growth_rates_F,
                individuals_per_zone_F,
                individuals_per_zone_M,
                5,
                6,
                1,
                10,
                active_F,
                active_M,
                inactive_F,
                inactive_M,
                buffer_M,
                inactive_locations_F,
                inactive_locations_M,
                inactive_genotypes_F,
                inactive_genotypes_M,
                inactive_mitochondria_F,
                inactive_mitochondria_M)
        end
    end

    # initializes new collection of population data variables by replacing the old individuals with their offspring and updating the growth rates
    # individuals in an inactive zone remain fixed and the boundaries for active/inactive zones are updated
    function PopulationData(pd,
        genotypes_daughters,
        genotypes_sons,
        mitochondria_daughters,
        mitochondria_sons,
        locations_daughters,
        locations_sons,
        expand_left,
        expand_right,
        generation,
        competition_trait_loci,
        K_total,
        sigma_comp,
        intrinsic_R,
        ecolDiff)

        new_active_F = [] # keeps track of inactive females that will be active in the next generation
        new_active_M = [] # keeps track of inactive males that will be active in the next generation
        buffer_M = pd.buffer_M
        individuals_per_zone_F = pd.individuals_per_zone_F
        individuals_per_zone_M = pd.individuals_per_zone_M
        left_boundary, right_boundary = pd.left_boundary, pd.right_boundary
        bottom_boundary, top_boundary = pd.bottom_boundary, pd.top_boundary
        inactive_locations_F, inactive_locations_M = pd.inactive_locations_F, pd.inactive_locations_M
        inactive_genotypes_F, inactive_genotypes_M = pd.inactive_genotypes_F, pd.inactive_genotypes_M
        inactive_mitochondria_F, inactive_mitochondria_M = pd.inactive_mitochondria_F, pd.inactive_mitochondria_M

        if generation % 10 == 0
            # every 10 generations the active and inactive zones are reset to reflect current distribution of genotypes

            # compile list of genotypes/locations of offspring and inactive individuals
            genotypes_F = vcat(genotypes_daughters, inactive_genotypes_F[pd.inactive_F])
            genotypes_M = vcat(genotypes_sons,
                inactive_genotypes_M[buffer_M[1]],
                inactive_genotypes_M[buffer_M[2]],
                inactive_genotypes_M[pd.inactive_M])

            mitochondria_F = vcat(mitochondria_daughters, inactive_mitochondria_F[pd.inactive_F])
            mitochondria_M = vcat(mitochondria_sons,
                inactive_mitochondria_M[buffer_M[1]],
                inactive_mitochondria_M[buffer_M[2]],
                inactive_mitochondria_M[pd.inactive_M])

            locations_F = vcat(locations_daughters, inactive_locations_F[pd.inactive_F])
            locations_M = vcat(locations_sons,
                inactive_locations_M[buffer_M[1]],
                inactive_locations_M[buffer_M[2]],
                inactive_locations_M[pd.inactive_M])

            # update inactive individuals to current population
            inactive_locations_F, inactive_locations_M = copy(locations_F), copy(locations_M)
            inactive_genotypes_F, inactive_genotypes_M = copy(genotypes_F), copy(genotypes_M)
            inactive_mitochondria_F, inactive_mitochondria_M = copy(mitochondria_F), copy(mitochondria_M)

            # group population into zones
            individuals_per_zone_F, individuals_per_zone_M = assign_zones(locations_F, locations_M)

            left_boundary, right_boundary, bottom_boundary, top_boundary = find_active_zone_boundaries(individuals_per_zone_F,
                individuals_per_zone_M,
                genotypes_F,
                genotypes_M,
                mitochondria_F,
                mitochondria_M)

            active_zones = [collect(left_boundary:right_boundary), collect(bottom_boundary:top_boundary)]

            dead_zones = [[(1:left_boundary - 1); (right_boundary+1):10], (1:10)]

            # compile list of new inactive individuals
            if dead_zones[1] != [] && dead_zones[2] !=[]
                inactive_F = reduce(vcat, individuals_per_zone_F[dead_zones[1], dead_zones[2]])
                inactive_M = reduce(vcat, individuals_per_zone_M[dead_zones[1], dead_zones[2]])
            else
                inactive_F = []
                inactive_M = []
            end

            # compile list of new active individuals
            if active_zones != []
                active_F = reduce(vcat, individuals_per_zone_F[active_zones[1], active_zones[2]])
                active_M = reduce(vcat, individuals_per_zone_M[active_zones[1], active_zones[2]])
            else
                active_F = []
                active_M = []
            end

            buffer_M = [[], []]
        else
            # expands active zone and updates buffer if necessary
            if expand_left && left_boundary > 1 # expands the left boundary if not already at edge of range
                left_boundary -= 1
                new_active_F = vcat(individuals_per_zone_F[left_boundary, :]...)
                new_active_M = vcat(individuals_per_zone_M[left_boundary, :]...)
                if left_boundary > 1
                    buffer_M[1] = vcat(individuals_per_zone_M[left_boundary-1, :]...)
                else
                    buffer_M[1] = []
                end
            end

            if expand_right && right_boundary < 10 # expands the right boundary if not already at edge of range
                right_boundary += 1
                new_active_F = vcat(new_active_F, vcat(individuals_per_zone_F[right_boundary, :]...))
                new_active_M = vcat(new_active_M, vcat(individuals_per_zone_M[right_boundary, :]...))
                if right_boundary < 10
                    buffer_M[2] = vcat(individuals_per_zone_M[right_boundary+1, :]...)
                else
                    buffer_M[2] = []
                end
            end

            new_active_M = vcat(new_active_M, buffer_M[1], buffer_M[2]) # add males in the buffer zone to the list of active males

            try
                # updates list of inactive individuals
                inactive_F = setdiff(pd.inactive_F, new_active_F)
                inactive_M = setdiff(pd.inactive_M, new_active_M)

                # creates new lists of genotypes and locations so that the active individuals are first and the inactive individuals last 
                genotypes_F = vcat(genotypes_daughters, inactive_genotypes_F[new_active_F], inactive_genotypes_F[inactive_F])
                genotypes_M = vcat(genotypes_sons, inactive_genotypes_M[new_active_M], inactive_genotypes_M[inactive_M])

                mitochondria_F = vcat(mitochondria_daughters, inactive_mitochondria_F[new_active_F], inactive_mitochondria_F[inactive_F])
                mitochondria_M = vcat(mitochondria_sons, inactive_mitochondria_M[new_active_M], inactive_mitochondria_M[inactive_M])

                locations_F = vcat(locations_daughters, inactive_locations_F[new_active_F], inactive_locations_F[inactive_F])
                locations_M = vcat(locations_sons, inactive_locations_M[new_active_M], inactive_locations_M[inactive_M])
            catch e
                println(eachindex(locations_F))
                println(eachindex(locations_M))

                @testset "find_error" begin
                    num_inactive_F = reduce(+, map(length, individuals_per_zone_F))
                    num_inactive_M = reduce(+, map(length, individuals_per_zone_M))

                    @test num_inactive_F == length(inactive_locations_F)
                    @test num_inactive_M == length(inactive_locations_M)
                    @test maximum(vcat(0, new_active_F)) <= length(inactive_locations_F)
                    @test maximum(vcat(0, new_active_M)) <= length(inactive_locations_M)
                    @test maximum(vcat(0, inactive_F)) <= length(inactive_locations_F)
                    @test maximum(vcat(0, inactive_M)) <= length(inactive_locations_M)
                end
            end

            # updates list of active individuals
            num_active_F = length(locations_daughters) + length(new_active_F)
            num_active_M = length(locations_sons) + length(new_active_M)

            active_F = collect(1:num_active_F)
            active_M = collect(1:num_active_M)
        end

        # calculate new competition traits
        competition_traits_F = calc_traits_additive(genotypes_F, competition_trait_loci)
        competition_traits_M = calc_traits_additive(genotypes_M, competition_trait_loci)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1

        ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)


        # calculate new growth rates
        growth_rates_F = calculate_growth_rates(ind_useResourceA_F,
            ind_useResourceB_F,
            ind_useResourceA_M,
            ind_useResourceB_M,
            locations_F,
            locations_M,
            K_total,
            sigma_comp,
            intrinsic_R,
            active_F)

        new(genotypes_F,
            genotypes_M,
            locations_F,
            locations_M,
            mitochondria_F,
            mitochondria_M,
            growth_rates_F,
            individuals_per_zone_F,
            individuals_per_zone_M,
            left_boundary,
            right_boundary,
            bottom_boundary,
            top_boundary,
            active_F,
            active_M,
            inactive_F,
            inactive_M,
            buffer_M,
            inactive_locations_F,
            inactive_locations_M,
            inactive_genotypes_F,
            inactive_genotypes_M,
            inactive_mitochondria_F,
            inactive_mitochondria_M)
    end

    # initializes new collection of population data variables by replacing the old individuals with their offspring and updating the growth rates
    function PopulationData(genotypes_daughters,
        genotypes_sons,
        mitochondria_daughters,
        mitochondria_sons,
        locations_daughters::Vector{Location},
        locations_sons::Vector{Location},
        competition_trait_loci,
        K_total,
        sigma_comp,
        intrinsic_R,
        ecolDiff)

        # calculate new competition traits
        competition_traits_F = calc_traits_additive(genotypes_daughters, competition_trait_loci)
        competition_traits_M = calc_traits_additive(genotypes_sons, competition_trait_loci)

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1

        ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M = calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)


        # calculate new growth rates
        growth_rates_F = calculate_growth_rates(ind_useResourceA_F,
            ind_useResourceB_F,
            ind_useResourceA_M,
            ind_useResourceB_M,
            locations_daughters,
            locations_sons,
            K_total,
            sigma_comp,
            intrinsic_R,
            1:length(locations_daughters))

        new(genotypes_daughters,
            genotypes_sons,
            locations_daughters,
            locations_sons,
            mitochondria_daughters,
            mitochondria_sons,
            growth_rates_F)
    end
end

# updates the genotypes, locations, and growth rates with the values for the next generation
function update_population(pd,
    genotypes_daughters,
    genotypes_sons,
    mitochondria_daughters,
    mitochondria_sons,
    locations_daughters::Vector{Location},
    locations_sons::Vector{Location},
    expand_left,
    expand_right,
    generation,
    competition_trait_loci,
    K_total,
    sigma_comp,
    intrinsic_R,
    ecolDiff,
    optimize)

    if optimize
        # updates the population with active/inactive zones
        return PopulationData(pd,
            genotypes_daughters,
            genotypes_sons,
            mitochondria_daughters,
            mitochondria_sons,
            locations_daughters,
            locations_sons,
            expand_left,
            expand_right,
            generation,
            competition_trait_loci,
            K_total,
            sigma_comp,
            intrinsic_R,
            ecolDiff)
    else
        # updates the population without active/inactive zones
        return PopulationData(genotypes_daughters,
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
    end
end

# Generate breeding locations of individuals
function generate_initial_locations(pop0_starting_N_half, pop1_starting_N_half, starting_range_pop0, starting_range_pop1)
    locations_F_pop0 = Vector{Location}(undef, pop0_starting_N_half)
    locations_M_pop0 = Vector{Location}(undef, pop0_starting_N_half)
    locations_F_pop1 = Vector{Location}(undef, pop1_starting_N_half)
    locations_M_pop1 = Vector{Location}(undef, pop1_starting_N_half)

    for i in 1:pop0_starting_N_half
        locations_F_pop0[i] = Location(starting_range_pop0)
        locations_M_pop0[i] = Location(starting_range_pop0)
    end

    for i in 1:pop1_starting_N_half
        locations_F_pop1[i] = Location(starting_range_pop1)
        locations_M_pop1[i] = Location(starting_range_pop1)
    end

    return [locations_F_pop0; locations_F_pop1], [locations_M_pop0; locations_M_pop1]
end

# arranges the indices of each individual by location
# zone 1 is [0,1), zone 2 is [1, 2), etc.
# list of indices for a zone is stored at individuals_per_zone_F[zone_x_coordinate][zone_y_coordinate]
function assign_zones(locations_F, locations_M)
    individuals_per_zone_F = Matrix{Vector{Integer}}(undef, 10, 10)
    individuals_per_zone_M = Matrix{Vector{Integer}}(undef, 10, 10)

    # creates an empty vector to store the list of indices at each entry in the 10x10 matrix
    for i in 1:10
        for j in 1:10
            individuals_per_zone_F[i, j] = Vector{Integer}(undef, 0)
            individuals_per_zone_M[i, j] = Vector{Integer}(undef, 0)
        end
    end

    # assigns females
    for indv in eachindex(locations_F)
        zone_x = trunc(Int, 10 * locations_F[indv].x) + 1
        zone_y = trunc(Int, 10 * locations_F[indv].y) + 1
        push!(individuals_per_zone_F[zone_x, zone_y], indv)
    end

    #assigns males
    for indv in eachindex(locations_M)
        zone_x = trunc(Int, 10 * locations_M[indv].x) + 1
        zone_y = trunc(Int, 10 * locations_M[indv].y) + 1
        push!(individuals_per_zone_M[zone_x, zone_y], indv)
    end

    return individuals_per_zone_F, individuals_per_zone_M
end

# finds which parts of the range contain only individuals of one genotype
function find_active_zone_boundaries(indv_per_zone_F, indv_per_zone_M, genotypes_F, genotypes_M, mitochondria_F, mitochondria_M)
    dead_zones = []
    left_boundary, bottom_boundary = 11,11
    right_boundary, top_boundary = 0,0

    for x_zone in 1:10
        females = vcat(indv_per_zone_F[x_zone, :]...)
        males = vcat(indv_per_zone_M[x_zone, :]...)

        genotypes = [genotypes_F[females]; genotypes_M[males]]
        mitochondria = [mitochondria_F[females]; mitochondria_M[males]]

        if length(mitochondria) == 0 || maximum(map(maximum, genotypes)) > 0 || maximum(mitochondria) > 0
            left_boundary = max(x_zone, 1)
            break
        end
    end
    for x_zone in 10:-1:left_boundary
        females = vcat(indv_per_zone_F[x_zone, :]...)
        males = vcat(indv_per_zone_M[x_zone, :]...)
        
        genotypes = [genotypes_F[females]; genotypes_M[males]]
        mitochondria = [mitochondria_F[females]; mitochondria_M[males]]

        if length(mitochondria) == 0 || minimum(map(minimum, genotypes)) < 1 || minimum(mitochondria) < 1
            right_boundary = min(x_zone, 10)
            break
        end
    end

    print(string("left: ", left_boundary))
    print(string("right: ", right_boundary))
    return left_boundary, right_boundary, 1, 10
end

# calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
function calculate_ind_useResource(competition_traits_F, competition_traits_M, ecolDiff)
    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    competAbility = (1 - ecolDiff) / 2    # equals 0 when ecolDiff = 1 

    ind_useResourceA_F = competAbility .+ ((1 .- competition_traits_F) .* ecolDiff)
    ind_useResourceB_F = competAbility .+ (competition_traits_F .* ecolDiff)
    ind_useResourceA_M = competAbility .+ ((1 .- competition_traits_M) .* ecolDiff)
    ind_useResourceB_M = competAbility .+ (competition_traits_M .* ecolDiff)

    return ind_useResourceA_F, ind_useResourceB_F, ind_useResourceA_M, ind_useResourceB_M
end

# calculates the ideal density (assuming normal distribution) at the location of each female
function get_ideal_densities(K_total, sigma_comp, locations_F)
    function erf(x) # approximation of the error function
        return tanh((sqrt(pi) * log(2)) * x)
    end

    function calc_ideal_density(location) # integral from 0 to 1 of K_total*exp(-(x-focal_location)^2/(2*sigma_comp^2)) with respect to x
        x_density = K_total * sqrt(pi / 2) * sigma_comp * (erf((1 - location.x) / (sqrt(2) * sigma_comp)) + erf(location.x / (sqrt(2) * sigma_comp)))
        y_density = K_total * sqrt(pi / 2) * sigma_comp * (erf((1 - location.y) / (sqrt(2) * sigma_comp)) + erf(location.y / (sqrt(2) * sigma_comp)))

        return x_density + y_density
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
function calculate_growth_rates(ind_useResourceA_F,
    ind_useResourceB_F,
    ind_useResourceA_M,
    ind_useResourceB_M,
    locations_F::Vector{Location},
    locations_M,
    K_total,
    sigma_comp,
    intrinsic_R,
    active_F)
    # in spatial model, calculate growth rates based on local resource use

    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    ideal_densities_at_locations_F = get_ideal_densities(K_total, sigma_comp, locations_F[active_F]) # this applies the above function to each geographic location

    ideal_densities_at_locations_F_resourceA = ideal_densities_at_locations_F ./ 2
    ideal_densities_at_locations_F_resourceB = ideal_densities_at_locations_F_resourceA

    # assume both resources have same constant density across range

    # determine local resource use for each location across range
    ind_useResourceA_all = [ind_useResourceA_F; ind_useResourceA_M]
    ind_useResourceB_all = [ind_useResourceB_F; ind_useResourceB_M]
    ind_locations_real = [locations_F; locations_M]

    function get_useResourceA_density_real(focal_location) # this function calculates local density according to a normal curve, weighted by individual resource use
        sum(ind_useResourceA_all .* exp.(-get_squared_distances(ind_locations_real, focal_location) ./ (2 * (sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end

    real_densities_at_locations_F_resourceA = map(get_useResourceA_density_real, locations_F[active_F]) # this applies the above function to each geographic location

    function get_useResourceB_density_real(focal_location) # do the same for resource B
        sum(ind_useResourceB_all .* exp.(-get_squared_distances(ind_locations_real, focal_location) ./ (2 * (sigma_comp^2))))
    end
    real_densities_at_locations_F_resourceB = map(get_useResourceB_density_real, locations_F[active_F]) # this applies the above function to each geographic location 

    # calculate local growth rates due to each resource (according to discrete time logistic growth equation)

    local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_locations_F_resourceA ./ (ideal_densities_at_locations_F_resourceA .+ ((real_densities_at_locations_F_resourceA) .* (intrinsic_R - 1)))
    local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_locations_F_resourceB ./ (ideal_densities_at_locations_F_resourceB .+ ((real_densities_at_locations_F_resourceB) .* (intrinsic_R - 1)))

    growth_rateA = ind_useResourceA_F[active_F] .* local_growth_rates_resourceA
    growth_rateB = ind_useResourceB_F[active_F] .* local_growth_rates_resourceB

    growth_rates = growth_rateA .+ growth_rateB

    return Dict(active_F .=> growth_rates)
end

# This function sets up the genotypes of the starting population
# in a 3D array, where:
# rows (D1) are alleles (row 1 from mother, row 2 from father),
# columns (D2) are loci, 
# pages (D3) are individuals.  
function generate_genotype_array(N_pop0::Integer, N_pop1::Integer, loci::Integer)
    total_N = N_pop0 + N_pop1
    genotypes = vcat([fill(0, 2, loci) for i in 1:N_pop0], [fill(1, 2, loci) for j in (N_pop0+1):total_N])

    return genotypes
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

# This function determines breeding location of one individual,
# based on a normal distribution with width sigma_disp,
# centred on birth location. Constrained to be within range. 
# geographic_limits should be a vector with two numbers.
function disperse_individual(female_location, sigma_disp::Real, geographic_limits::Vector{})::Float32
    while true
        new_location = female_location + (sigma_disp * randn())
        if (new_location >= geographic_limits[1]) & (new_location <= geographic_limits[2]) # checks if location is in range
            return new_location
        end
    end
end

# creates a plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function plot_population(pd, optimize, functional_loci_range)
    if optimize
        genotypes_active = [pd.genotypes_F[pd.active_F]; pd.genotypes_M[pd.active_M]]
        locations_active = [pd.locations_F[pd.active_F]; pd.locations_M[pd.active_M]]
        genotypes_inactive = [pd.inactive_genotypes_F[pd.inactive_F]; pd.inactive_genotypes_M[pd.inactive_M]]
        locations_inactive = [pd.inactive_locations_F[pd.inactive_F]; pd.inactive_locations_M[pd.inactive_M]]
        mitochondria_active = [pd.mitochondria_F[pd.active_F]; pd.mitochondria_M[pd.active_M]]
        mitochondria_inactive = [pd.inactive_mitochondria_F[pd.inactive_F]; pd.inactive_mitochondria_M[pd.inactive_M]]

        create_new_plot(calc_traits_additive(genotypes_active, functional_loci_range),
            mitochondria_active,
            locations_active,
            calc_traits_additive(genotypes_inactive, functional_loci_range),
            mitochondria_inactive,
            locations_inactive)
    else
        create_new_plot(calc_traits_additive([pd.genotypes_F; pd.genotypes_M], functional_loci_range),
            [pd.mitochondria_F; pd.mitochondria_M],
            [pd.locations_F; pd.locations_M])
    end
end

# updates the plot of locations and hybrid indices (plotting handled by the Plot_Data module)
function update_plot(pd, generation, optimize, functional_loci_range)
    if optimize 
        genotypes_active = [pd.genotypes_F[pd.active_F]; pd.genotypes_M[pd.active_M]]
        locations_active = [pd.locations_F[pd.active_F]; pd.locations_M[pd.active_M]]
        genotypes_inactive = [pd.inactive_genotypes_F[pd.inactive_F]; pd.inactive_genotypes_M[pd.inactive_M]]
        locations_inactive = [pd.inactive_locations_F[pd.inactive_F]; pd.inactive_locations_M[pd.inactive_M]]
        mitochondria_active = [pd.mitochondria_F[pd.active_F]; pd.mitochondria_M[pd.active_M]]
        mitochondria_inactive = [pd.inactive_mitochondria_F[pd.inactive_F]; pd.inactive_mitochondria_M[pd.inactive_M]]

        update_population_plot(calc_traits_additive(genotypes_active, functional_loci_range),
            mitochondria_active,
            locations_active,
            calc_traits_additive(genotypes_inactive, functional_loci_range),
            mitochondria_inactive,
            locations_inactive,
            generation)
    else
        genotypes = [pd.genotypes_F; pd.genotypes_M]
        locations = [pd.locations_F; pd.locations_M]
        mitochondria = [pd.mitochondria_F; pd.mitochondria_M]
        update_population_plot(calc_traits_additive(genotypes, functional_loci_range),
            mitochondria,
            locations,
            generation)
    end
end
end