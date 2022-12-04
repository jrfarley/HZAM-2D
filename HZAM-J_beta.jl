#= HZAM-J_beta.jl
Hybrid Zone with Assortative Mating, Julia implementation,
by Darren Irwin.
Starting this file on 30Jan2022 to develop an integrated HZAM model that 
combines features of HZAM_Sym_Julia_V1.jl and HZAM_release_2.0.R into 

If you use this code, please cite the paper below, and possibly my GitHub or Dryad repository 
where you obtained this code:

Irwin, D., and Schluter, D. Hybridization and the coexistence of species. bioRxiv 2021.04.04.438369; doi: https://doi.org/10.1101/2021.04.04.438369 
=#

#= You need to make sure these packages below are added to your Julia environment.
To add e.g. the package Distributions, 
type "]" to get the "pkg>" prompt, then type e.g. "add Distributions";
or run the commands below: 
=#
# import Pkg; 
# Pkg.add("Distributions") 
# Pkg.add("Statistics") 
# Pkg.add("JLD2")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")
# Pkg.add("CategoricalArrays")
# Pkg.add("Colors")
# Pkg.add("ColorSchemes")
# Pkg.add("CategoricalArrays")
# Pkg.add("Makie") 
# Pkg.add("CairoMakie")
# Pkg.add("GLMakie")
# Pkg.add("LsqFit") 

using Distributions # needed for "Poisson" function
using Statistics  # needed for "mean" function
using JLD2 # needed for saving / loading data in Julia format
using CSV # for saving in csv format
using DataFrames # for converting data to save as csv
using LsqFit

# for plotting:
# using Plots
# gr()  # use GR backend for graphs
using CategoricalArrays
using Colors, ColorSchemes
import ColorSchemes.plasma
# using Plots.PlotMeasures  # needed for plot margin adjustment
using GLMakie

# to start Julia with multiple threads, type in terminal e.g.:
# julia --threads 4
# To check, type in Julia: Threads.nthreads()
# 4

# set up functions (should not define functions repeatedly in loop, as causes re-compilation, slows things)

# This function sets up the genotypes of the starting population
# in a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
# columns (D2) are loci, and pages (D3) are individuals
function generate_genotype_array(N_pop0,N_pop1,loci)::Array{Int8, 3}
    total_N = N_pop0 + N_pop1  
    genotypes = Array{Int8, 3}(undef, 2, loci, total_N) # The "Int8" is the type (8-bit integer), and "undef" means an unitialized array, so values are meaningless
    genotypes[:,:,1:N_pop0] .= 0  # assigns genotypes of pop01
    genotypes[:,:,(N_pop0+1):total_N] .= 1  # assigns genotypes of pop1
    return genotypes
end

function calc_traits_additive(genotypes::Array{Int, 3})::Vector{Float32}
    N = size(genotypes, 3) 
    traits = Vector{Float32}(undef, N) # Float32 should be enough precision; memory saving compared to Float64
    for i in 1:N
        traits[i] = mean(genotypes[:,:,i])
    end
    return traits
end

function get_survival_fitnesses_epistasis(genotypes, w_hyb, beta=1)
    survival_HI = calc_traits_additive(genotypes)
    epistasis_fitnesses = 1 .- (1 - w_hyb) .* (4 .* survival_HI .* (1 .- survival_HI)).^beta
    return epistasis_fitnesses 
end

function get_survival_fitnesses_hetdisadvantage(genotypes, w_hyb)
    N = size(genotypes, 3)
    num_loci = size(genotypes, 2)
    s_per_locus = 1 - w_hyb ^ (1/num_loci)  # loss in fitness due to each heterozygous locus 
    num_hetloci = Array{Int16, 1}(undef, N)
    for ind in 1:N  # count number of het loci per individual
        num_hetloci[ind] = sum(genotypes[1,:,ind] .!= genotypes[2,:,ind])
    end
    hetdisadvantage_fitnesses = (1 - s_per_locus) .^ num_hetloci
    return hetdisadvantage_fitnesses
end

function disperse_individual(start_location, sigma_disp, geographic_limits)::Float32
    while true
        new_location = start_location + (sigma_disp * randn())
        if (new_location >= geographic_limits[1]) & (new_location <= geographic_limits[2]) # checks if location is in range
            return new_location 
        end
    end 
end

# These next two functions define the way potential male mates are chosen
# When in sympatric model, choose random male (note the second and third arguments not used but allows function to be called in same way as in spatial model):
function choose_random_male(elig_M, locations_M, focal_location)
    focal_male = splice!(elig_M, rand(eachindex(elig_M)))
    return focal_male, elig_M
end
# When in spatial model, choose closest male:
function choose_closest_male(elig_M, locations_M, focal_location)
    focal_male = splice!(elig_M, argmin(abs.(locations_M[elig_M] .- focal_location)))
    return focal_male, elig_M
end

function sigmoid(x, p)  # where p[1] is cline centre, and p[2] is maxSlope 
    1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

# for adding jitter (small random shifts in position, to better visualize overlapping points)
function jitter(n::Vector{Float32}, factor=0.02) 
    n .+ (0.5 .- rand(length(n))) .* factor
end

function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int = 1000, max_generations::Int = 1000, 
    total_loci::Int = 6, female_mating_trait_loci = 1:3, male_mating_trait_loci = 1:3,
    competition_trait_loci = 1:3, hybrid_survival_loci = 1:3, neutral_loci = 4:6,
    survival_fitness_method::String = "epistasis", per_reject_cost = 0,
    starting_pop_ratio = 1.0, sympatry = false, geographic_limits = [0, 1], 
    starting_range_pop0 = [0.0, 0.48], starting_range_pop1 = [0.52, 1.0],
    sigma_disp = 0.01, sigma_comp = 0.01,
    do_plot = true, plot_int = 10)

    # specify ecological resource competitive abilities for two resources A and B 
    # ecolDiff = 1.0 # this is "E" in the paper 
    competAbility_useResourceA_pop0 = (1 + ecolDiff)/2    # equals 1 when ecolDiff = 1   
    competAbility_useResourceB_pop0 = 1 - competAbility_useResourceA_pop0
    competAbility_useResourceA_pop1 = (1 - ecolDiff)/2   # equals 0 when ecolDiff = 1
    competAbility_useResourceB_pop1 = 1 - competAbility_useResourceA_pop1

    # set up carying capacities on each resource, and starting pop sizes of each species
    K_A = K_total / 2  # EVEN NUMBER; carrying capacity (on resource alpha) of entire range (for two sexes combined), regardless of species 
    K_B = K_total / 2   # EVEN NUMBER; carrying capacity (on resource beta) of entire range (for two sexes combined), regardless of species

    if sympatry  # if sympatry is true, then set both pop sizes according to full range 
        starting_range_pop0 = [0.0, 1.0] 
        starting_range_pop1 = [0.0, 1.0]
        pop0_starting_N = K_A   # starting N of species 0
        pop0_starting_N_half = Int(pop0_starting_N/2)  # The "Int" is to ensure no decimal
        pop1_starting_N = Int(round(starting_pop_ratio * K_B))   # starting N of species 1, which can be lower if starting_pop_ratio is below 1)
        pop1_starting_N_half = Int(pop1_starting_N/2)
    else  # sympatry is false, then set pop sizes assuming no range overlap--this is why the "2" is in the formulae below  
        pop0_starting_N = (2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop0) + (K_B * competAbility_useResourceB_pop0)) * (starting_range_pop0[2] - starting_range_pop0[1]) / (geographic_limits[2] - geographic_limits[1]) 
        pop0_starting_N_half = Int(pop0_starting_N/2)  # starting number for each sex
        pop1_starting_N = (2 - ecolDiff) * ((K_A * competAbility_useResourceA_pop1) + (K_B * competAbility_useResourceB_pop1)) * (starting_range_pop1[2] - starting_range_pop1[1]) / (geographic_limits[2] - geographic_limits[1])
        pop1_starting_N_half = Int(pop1_starting_N/2)
    end
    
    # WHEN ECOLDIFF = 1, THEN SHOULD NOT HAVE THE 2 (1 INSTEAD)
    # WHEN ECOLDIFF = 0, THEN SHOULD HAVE THE 2
    # SO 2-ecolDiff ?

    # get the chosen survival fitness function
    if survival_fitness_method == "epistasis"
        get_survival_fitnesses = get_survival_fitnesses_epistasis
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        get_survival_fitnesses = get_survival_fitnesses_hetdisadvantage
        short_survFitnessMethod = "Het"
    else
        println("ERROR--no survival fitness method chosen--should be either epistasis or hetdisadvantage")
    end

    # add up and check the total loci
    total_functional_loci = max(maximum(female_mating_trait_loci), maximum(male_mating_trait_loci), maximum(competition_trait_loci), maximum(hybrid_survival_loci))
    functional_loci_range = 1:total_functional_loci
    num_neutral_loci = length(Vector(neutral_loci))
    if total_functional_loci + num_neutral_loci ≠ total_loci
        println("#### WARNING: Please examine your loci numbers and indices, as they don't all match up ####")
    end
    
    # convert S_AM to pref_ratio (for use in math below)
    if S_AM == 1
        S_AM_for_math = 1 + 10^(-15)
    elseif S_AM == Inf
        S_AM_for_math = 10^(15)
    else
        S_AM_for_math = S_AM
    end
    pref_SD = sqrt(-1 / (2 * log(1 / S_AM_for_math)))  # width of female acceptance curve for male trait

    # Generate genotype array for population of females:
    # this is a 3D array, where rows (D1) are alleles (row 1 from mother, row 2 from father),
    # columns (D2) are loci, and pages (D3) are individuals
    genotypes_F = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    genotypes_M = generate_genotype_array(pop0_starting_N_half, pop1_starting_N_half, total_loci)
    # functional loci are first, followed by neutral loci

    # Generate breeding locations of individuals (vector in same order of individuals as D3 of the genotype array)
    locations_F_pop0 = Array{Float32, 1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1]) 
    locations_F_pop1 = Array{Float32, 1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_F = [locations_F_pop0; locations_F_pop1]
    locations_M_pop0 = Array{Float32, 1}((rand(pop0_starting_N_half) .* (starting_range_pop0[2] - starting_range_pop0[1])) .+ starting_range_pop0[1]) 
    locations_M_pop1 = Array{Float32, 1}((rand(pop1_starting_N_half) .* (starting_range_pop1[2] - starting_range_pop1[1])) .+ starting_range_pop1[1])
    locations_M = [locations_M_pop0; locations_M_pop1]
    
    # set up expected local densities, based on geographically even distribution of individuals at carrying capacity
    spaced_locations = collect(Float32, geographic_limits[1]:0.001:geographic_limits[2])
    ind_locations_if_even_at_K = range(geographic_limits[1], geographic_limits[2], length=K_total)
    function get_density_if_even_at_K(focal_location) # this function calculates local density according to a normal curve
        sum(exp.(-((ind_locations_if_even_at_K .- focal_location).^2)./(2*(sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
    end
    ideal_densities_at_spaced_locations = map(get_density_if_even_at_K, spaced_locations) # this applies the above function to each geographic location
    # assume both resources have same constant density across range
    ideal_densities_at_spaced_locations_resourceA = ideal_densities_at_spaced_locations ./ 2
    ideal_densities_at_spaced_locations_resourceB = ideal_densities_at_spaced_locations ./ 2 

    #plot(spaced_locations, ideal_densities_at_spaced_locations_resourceA)

    extinction = false  # if extinction happens later this will be set true

    # Pick function for choosing potential male mate (these two functions defined above) 
    if sympatry # if sympatry, then choose random male
        pick_potential_mate = choose_random_male
    else # if space matters, then choose closest male
        pick_potential_mate = choose_closest_male
    end
    
    if do_plot
        functionalLoci_HI_all_inds = [calc_traits_additive(genotypes_F[:, functional_loci_range, :]); calc_traits_additive(genotypes_M[:, functional_loci_range, :])]
        fontsize_theme = Theme(fontsize = 60); set_theme!(fontsize_theme)  # this sets the standard font size
        fig = Figure(resolution=(1800, 1200), figure_padding = 60)
        ax = Axis(fig[1, 1], xlabel = "location", ylabel = "hybrid index", title = string("HZAM simulation, generation = ", 0), xticklabelsize = 45, yticklabelsize = 45, titlegap = 30)
        xlims!(-0.03, 1.03)
        ylims!(-0.03, 1.03)
        points = scatter!(ax, [locations_F; locations_M], jitter(functionalLoci_HI_all_inds), color = (:blue, 0.5))
        initial_par = [0., 1.]  # next line will start search for centre and slope with these values
        fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
        sigmoid_line = lines!(ax, spaced_locations, sigmoid(spaced_locations, fit.param), color = (:blue, 0.25), linewidth = 20)
        display(fig)
    end

    # loop throught the generations
    for generation in 1:max_generations

        # Prepare for mating and reproduction
        N_F = size(genotypes_F, 3)
        N_M = size(genotypes_M, 3)

        println("generation: ", generation, "; individuals: ", N_F + N_M)
 
        # calculate mating trait values (T) from genotypes
        female_mating_traits = calc_traits_additive(genotypes_F[:, female_mating_trait_loci, :])
        male_mating_traits = calc_traits_additive(genotypes_M[:, male_mating_trait_loci, :])

        # calculate ecological competition trait values from genotypes
        competition_traits_F = calc_traits_additive(genotypes_F[:, competition_trait_loci, :])
        competition_traits_M = calc_traits_additive(genotypes_M[:, competition_trait_loci, :])

        # calculate individual contributions to resource use, according to linear gradient between use of species 0 and species 1
        ind_useResourceA_F = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_F) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
        ind_useResourceB_F = competAbility_useResourceB_pop0 .+ (competition_traits_F .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))
        ind_useResourceA_M = competAbility_useResourceA_pop1 .+ ((1 .- competition_traits_M) .* (competAbility_useResourceA_pop0 - competAbility_useResourceA_pop1))
        ind_useResourceB_M = competAbility_useResourceB_pop0 .+ (competition_traits_M .* (competAbility_useResourceB_pop1 - competAbility_useResourceB_pop0))
         
        if sympatry
            # sum up the global resource use over all individuals
            total_useResourceA = sum(ind_useResourceA_F) + sum(ind_useResourceA_M)
            total_useResourceB = sum(ind_useResourceB_F) + sum(ind_useResourceB_M)
            # calculate global growth rates due to each resource (according to discrete time logistic growth equation)
            growth_rate_resourceA = (intrinsic_R * K_A) / (K_A + ((total_useResourceA) * (intrinsic_R - 1)))
            growth_rate_resourceB = (intrinsic_R * K_B) / (K_B + ((total_useResourceB) * (intrinsic_R - 1)))
        else  # in spatial model, calculate growth rates based on local resource use
            # determine local resource use for each location across range
            ind_locations_real = [locations_F; locations_M]
            ind_useResourceA_all = [ind_useResourceA_F; ind_useResourceA_M]
            ind_useResourceB_all = [ind_useResourceB_F; ind_useResourceB_M]  
            function get_useResourceA_density_real(focal_location) # this function calculates local density according to a normal curve, weighted by individual resource use
                sum(ind_useResourceA_all .* exp.(-((ind_locations_real .- focal_location).^2)./(2*(sigma_comp^2)))) # because this function is within a function, it can use the variables within the larger function in its definition
            end
            real_densities_at_spaced_locations_resourceA = map(get_useResourceA_density_real, spaced_locations) # this applies the above function to each geographic location
            function get_useResourceB_density_real(focal_location) # do the same for resource B
                sum(ind_useResourceB_all .* exp.(-((ind_locations_real .- focal_location).^2)./(2*(sigma_comp^2))))
            end
            real_densities_at_spaced_locations_resourceB = map(get_useResourceB_density_real, spaced_locations) # this applies the above function to each geographic location 
            # calculate local growth rates due to each resource (according to discrete time logistic growth equation)
            local_growth_rates_resourceA = intrinsic_R .* ideal_densities_at_spaced_locations_resourceA ./ (ideal_densities_at_spaced_locations_resourceA .+ ((real_densities_at_spaced_locations_resourceA) .* (intrinsic_R - 1)))
            local_growth_rates_resourceB = intrinsic_R .* ideal_densities_at_spaced_locations_resourceB ./ (ideal_densities_at_spaced_locations_resourceB .+ ((real_densities_at_spaced_locations_resourceB) .* (intrinsic_R - 1)))
        end
        
        #plot(spaced_locations, local_growth_rates_resourceA)
        #plot(spaced_locations, local_growth_rates_resourceB)  
        
        #plot(spaced_locations, local_growth_rates)  
        
        ### NEED TO CAREFULLY PROOF THE ABOVE, PARTICULARLY IF ECOLDIF > 0


        # Set up structure to record number of matings per male (and female, although almost always 1 for females), 
        # to determine sexual selection due to HI class:
        matings_per_male = zeros(Int8, N_M)
        matings_per_female = zeros(Int8, N_F)

        # create empty data structures for keeping track of numbers offspring of parents
        daughters_per_mother = zeros(Int16, N_F)
        sons_per_mother = zeros(Int16, N_F)
        daughters_per_father = zeros(Int16, N_M)
        sons_per_father = zeros(Int16, N_M)

        # make empty arrays for storing genotypes of daughters and sons
        genotypes_daughters = Array{Int8,3}(undef, 2, total_loci, 0)
        genotypes_sons = Array{Int8,3}(undef, 2, total_loci, 0)
        
        # make empty arrays for storing locations of daughters and sons
        locations_daughters = Array{Float32, 1}(undef, 0)
        locations_sons = Array{Float32, 1}(undef, 0) 

        # create structures for recording indices of mother and father (for error checking, and potentially for tracking genealogies)
        # first column for mother index (3rd dim of genotypes_F) and second column for father index (3rd dim of genotypes_M) 
        daughter_parent_IDs = Array{Int}(undef, 0, 2)
        son_parent_IDs = Array{Int}(undef, 0, 2)

        # loop through mothers, mating and reproducing
        for mother in 1:N_F
            # initialize tracking variables
            mate = false # becomes true when female is paired with male
            rejects = 0 # will track number of rejected males (in cases there is cost--which there isn't in main HZAM-Sym paper)
            father = [] # will contain the index of the male mate
            # make vector of indices of eligible males
            elig_M = Vector(1:N_M)
            # determine male mate of female
            while mate == false
                # present female with random male (sympatric case) or closest male (spatial case), and remove him from list:
                focal_male, elig_M = pick_potential_mate(elig_M, locations_M, locations_F[mother])
                # compare male trait with female's trait (preference), and determine
                # whether she accepts; note that match_strength is determined by a
                # Gaussian, with a maximum of 1 and minimum of zero.
                match_strength = (exp(1)^((-(male_mating_traits[focal_male] - female_mating_traits[mother])^2) / (2 * (pref_SD^2))))
                if rand() < match_strength
                    # she accepts male, and mates
                    father = focal_male
                    matings_per_male[focal_male] += 1 # this adds 1 to the matings for that male
                    matings_per_female[mother] += 1
                    mate = true
                else
                    # she rejects male
                    rejects += 1
                    if length(elig_M) == 0
                        break
                    end
                end
            end
            # determine fitness cost due to mate search (number of rejected males)
            search_fitness = (1 - per_reject_cost)^rejects    # (in most of HZAM-Sym paper, per_reject_cost = 0)

            # determine fitness due to female use of available resources
            if sympatry
                growth_rate_of_focal_female = (ind_useResourceA_F[mother] * growth_rate_resourceA) + (ind_useResourceB_F[mother] * growth_rate_resourceB)
            else
                location_index = argmin(abs.(spaced_locations .- locations_F[mother]))
                local_growth_A = local_growth_rates_resourceA[location_index]
                local_growth_B = local_growth_rates_resourceB[location_index] 
                growth_rate_of_focal_female = (ind_useResourceA_F[mother] * local_growth_A) + (ind_useResourceB_F[mother] * local_growth_B)
            end

            #combine for total fitness:   
            reproductive_fitness = 2 * growth_rate_of_focal_female * search_fitness  # the 2 is because only females, not males, produce offspring
            # now draw the number of offspring from a poisson distribution with a mean of reproductive_fitness
            if isempty(father)
                offspring = 0  # if no mate (because all males were rejected), then no offspring
            else
                offspring = rand(Poisson(reproductive_fitness))
            end
            # if offspring, generate their genotypes and sexes
            if offspring >= 1
                for kid in 1:offspring
                    kid_info = Array{Int8,2}(undef, 2, total_loci) # place to store genotype of one offspring
                    # generate genotypes; for each locus (column), first row for allele from mother, second row for allele from father
                    for locus in 1:total_loci
                        kid_info[1, locus] = genotypes_F[rand([1 2]), locus, mother]  # for this locus, pick a random allele from the mother
                        kid_info[2, locus] = genotypes_M[rand([1 2]), locus, father]  # and from the father
                    end
                    # determine sex and location of kid
                    if rand() > 0.5 # kid is daughter
                        genotypes_daughters = cat(genotypes_daughters, kid_info, dims=3)
                        daughter_parent_IDs = cat(daughter_parent_IDs, [mother father], dims=1)
                        daughters_per_mother[mother] += 1
                        daughters_per_father[father] += 1
                        locations_daughters = [locations_daughters; disperse_individual(locations_F[mother], sigma_disp, geographic_limits)]
                    else # kid is son
                        genotypes_sons = cat(genotypes_sons, kid_info, dims=3)
                        son_parent_IDs = cat(son_parent_IDs, [mother father], dims=1)
                        sons_per_mother[mother] += 1
                        sons_per_father[father] += 1
                        locations_sons = [locations_sons; disperse_individual(locations_F[mother], sigma_disp, geographic_limits)]
                    end
                end
            end
        end # of loop through mothers

        # check if either no daughters or no sons, and end the simulation if so
        if (size(genotypes_daughters, 3) == 0) || (size(genotypes_sons, 3) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end

        # For someday: add in here the option of tracking fitness?

        # determine survival fitnesses of daughters due to epistasis

        survival_fitness_daughters = get_survival_fitnesses(genotypes_daughters[:, hybrid_survival_loci, :], w_hyb)
        daughters_survive = survival_fitness_daughters .> rand(length(survival_fitness_daughters))
        # same for sons:
        survival_fitness_sons = get_survival_fitnesses(genotypes_sons[:, hybrid_survival_loci, :], w_hyb)
        sons_survive = survival_fitness_sons .> rand(length(survival_fitness_sons))

        # check if either no surviving daughters or no surviving sons, and end the simulation if so
        if (sum(daughters_survive) == 0) || (sum(sons_survive) == 0)
            extinction = true # record an extinction of whole population (both "species")
            break # break out of current loop (this simulation) 
        end

        # assign surviving offspring to new adult population
        genotypes_F = genotypes_daughters[:, :, daughters_survive]
        locations_F = locations_daughters[daughters_survive] 
        genotypes_M = genotypes_sons[:, :, sons_survive]
        locations_M = locations_sons[sons_survive] 

        # update the plot
        if (do_plot && (generation % plot_int == 0))
            functionalLoci_HI_all_inds = [calc_traits_additive(genotypes_F[:, functional_loci_range, :]); calc_traits_additive(genotypes_M[:, functional_loci_range, :])] 
            # display(scatter([locations_F; locations_M], functionalLoci_HI_all_inds))
            # fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
            # lines!(spaced_locations, sigmoid(spaced_locations, fit.param))  # add the sigmoid fit to the plot
            delete!(ax, points)
            delete!(ax, sigmoid_line)
            points = scatter!([locations_F; locations_M], jitter(functionalLoci_HI_all_inds), color = (:blue, 0.5))
            fit = curve_fit(sigmoid, [locations_F; locations_M], functionalLoci_HI_all_inds, initial_par)
            sigmoid_line = lines!(spaced_locations, sigmoid(spaced_locations, fit.param), color = (:blue, 0.25), linewidth = 20)  # add the sigmoid fit to the plot
            ax.title = string("HZAM simulation, generation = ", generation) 
        end

    end # of loop through generations
    return genotypes_F, locations_F, genotypes_M, locations_M, extinction 
end
  

#### Run the actual simulation by calling the above function:

sim_results = run_one_HZAM_sim(0.9, 1000, 0, 1.1; # these values are 
                                # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total = 1000, max_generations = 1000,
    sigma_disp = 0.02, sympatry = false,
    sigma_comp = 0.1, do_plot = true, plot_int = 5)



#### The above is all that is needed for a single simulation.

####################################################################
# Below is additional things I have used for testing,
# and summarizing results from multiple simulations 




    #total_loci = 11,
    #female_mating_trait_loci = 1:5,
    #male_mating_trait_loci = 6:10,
    #competition_trait_loci = 6:10,
    #hybrid_survival_loci = 6:10,
    #neutral_loci = 11:11)








#------------------

    #total_loci = 20,
    #female_mating_trait_loci = 1:17,
    #male_mating_trait_loci = 1:17,
    #competition_trait_loci = 1:17,
    #hybrid_survival_loci = 1:17,
    #neutral_loci = 18:20)
functional_loci_range = 1:3
genotypes_F = sim_results[1]
genotypes_M = sim_results[2]
functional_HI_all_inds = [calc_traits_additive(genotypes_F[:, functional_loci_range, :]); calc_traits_additive(genotypes_M[:, functional_loci_range, :])]

# make a Makie plot

scatter(1:length(functional_HI_all_inds), functional_HI_all_inds)

# To test code inside the function:
w_hyb = 1
S_AM = 100000
ecolDiff = 1
intrinsic_R = 1.05
K_total = 5000 
max_generations = 200
total_loci = 6
female_mating_trait_loci = 1:3
male_mating_trait_loci = 1:3
competition_trait_loci = 1:3
hybrid_survival_loci = 1:3
neutral_loci = 4:6
survival_fitness_method = "epistasis"
per_reject_cost = 0
starting_pop_ratio = 1.0
sympatry = false
geographic_limits = [0, 1]
starting_range_pop0 = [0.0, 0.48]
starting_range_pop1 = [0.52, 1.0]
sigma_disp = 0.01
sigma_comp = 0.01
do_plot = true
plot_int = 1


# copy of the start of the function from above:
function run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;   # the semicolon makes the following optional keyword arguments  
    K_total::Int = 1000, max_generations::Int = 1000, 
    total_loci::Int = 6, female_mating_trait_loci = 1:3, male_mating_trait_loci = 1:3,
    competition_trait_loci = 1:3, hybrid_survival_loci = 1:3, neutral_loci = 4:6,
    survival_fitness_method::String = "epistasis", per_reject_cost = 0,
    starting_pop_ratio = 1.0, sympatry = false, geographic_limits = [0, 1], 
    starting_range_pop0 = [0.0, 0.48], starting_range_pop1 = [0.52, 1.0],
    sigma_disp = 0.01, sigma_comp = 0.01,
    do_plot = true, plot_int = 10)






#-------------------------
# functions used in Irwin & Schluter 2022 (the HZAM-sym model and graphing of results)

function run_HZAM_set(set_name::String, ecolDiff, intrinsic_R, replications;  # the semicolon makes the following optional keyword arguments 
    K_total::Int = 1000, max_generations::Int = 1000, 
    total_loci::Int = 6, female_mating_trait_loci = 1:3, male_mating_trait_loci = 1:3,
    competition_trait_loci = 1:3, hybrid_survival_loci = 1:3, neutral_loci = 4:6,
    survival_fitness_method::String = "epistasis", per_reject_cost = 0,
    starting_pop_ratio = 1.0)
    # ecolDiff should be from 0 to 1 (parameter "E" in the paper)
    # intrinsic_R is called "R" in the paper
    # replications should be somehting like "1:10" or just "1" for 1 replicate, or something like "2:5" to add replicates after 1 is done
    # K_total (default 1000) is the total carrying capacity, divided equally between K_alpha and K_beta
    # max_generations (default 1000) is the number of generations each simulation will run 
    # total_loci (default 6) is the total number of loci, regardless of their role (with indices 1:total_loci, referred to below)
    # Need to specify indices (columns) of four types of functional loci (can be the same). At least one should begin with index 1:
    # female_mating_trait_loci (default 1:3) is indices of the loci that determine the female mating trait
    # male_mating_trait_loci (default 1:3) is indices of the loci that determine the male mating trait
    # competition_trait_loci (default 1:3) is indices of the loci that determine the ecological trait (used in fitness related to resource use)
    # hybrid_survival_loci (default 1:3) is indices of the loci that determine survival probability of offspring to adulthood (can be viewed as incompatibilities and/or fitness valley based on ecology)
    # Specify indices (columns) of neutral loci (which have no effect on anything, just along for the ride):
    # neutral_loci (default 4:6) is indices of neutral loci (used for neutral measure of hybrid index; not used in the HZAM-sym paper)
    # per_reject_cost (default 0, can take values of 0 to 1) is fitness loss of female per male rejected (due to search time, etc.)
    # starting_pop_ratio (default 1) is the starting ratio of pop A to pop B

    save_outcomes_JL = true
    save_outcomes_csv = true  # whether to save the whole outcome array as csv files (with each rep as separate file)
    save_each_sim = false  # whether to save detailed data for each simulation

    # the set of hybrid fitnesses (w_hyb) values that will be run
    w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    # the set of assortative mating strengths (S_AM) that will be run
    S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]  # ratio of: probably of accepting homospecific vs. prob of accepting heterospecific

    if survival_fitness_method == "epistasis"
        short_survFitnessMethod = "Ep"
    elseif survival_fitness_method == "hetdisadvantage"
        short_survFitnessMethod = "Het"
    end

    # set up array of strings to record outcomes
    outcome_array = Array{String, 3}(undef, length(w_hyb_set), length(S_AM_set), length(replications))

    for k in 1:length(replications)  # loop through the replicate runs
        replicate_ID = replications[k]

        run_set_name = string(set_name,"_rep", replicate_ID)

        # Loop through the different simulation sets
        Threads.@threads for i in 1:length(w_hyb_set)
            for j in 1:length(S_AM_set)
                w_hyb = w_hyb_set[i]
                S_AM = S_AM_set[j]
                println("w_hyb = ", w_hyb, "; S_AM = ", S_AM)
            
                run_name = string("HZAM_animation_run", run_set_name, "_surv", short_survFitnessMethod, "_ecolDiff", ecolDiff, "_growthrate", intrinsic_R, "_K", K_total, "_FL", total_functional_loci, "_NL", num_neutral_loci, "_gen", max_generations, "_SC", per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM)
            
                # set up initial values for one simulation
                extinction = false
                outcome = []
                final_distribution = []
            
                # run one simulation by calling the function defined above:
                run_one_HZAM_sim(w_hyb, S_AM, ecolDiff, intrinsic_R;
                    K_total, max_generations,
                    total_loci, female_mating_trait_loci, male_mating_trait_loci,
                    competition_trait_loci, hybrid_survival_loci, neutral_loci,
                    survival_fitness_method, per_reject_cost,
                    starting_pop_ratio)
            
                # Record results of the one simulation
                functional_HI_all_inds = []
                species0_proportion = []
                species1_proportion = []
                HI_NL_all_inds = []
                species0_proportion_NL = []
                species1_proportion_NL = []
                if extinction  # whole simulation went extinct
                    outcome = "extinction"
                    if save_each_sim
                        @save string("simulation_data.", run_name, ".jld2") outcome
                    end
                else  # no complete extinction
                    # use trait loci to calculate HI of each individual
                    functional_HI_all_inds = [calc_traits_additive(genotypes_F[:, functional_loci_range, :]); calc_traits_additive(genotypes_M[:, functional_loci_range, :])]
                    # calculate proportion of all individuals who are species0 or species1 (defined as low and high 10% of HI distribution, respectively)
                    species0_proportion = sum(functional_HI_all_inds .< 0.1) / length(functional_HI_all_inds)
                    species1_proportion = sum(functional_HI_all_inds .> 0.9) / length(functional_HI_all_inds)
                    if species0_proportion >= 0.85 || species1_proportion >= 0.85
                        outcome = "one_species"
                    elseif (species0_proportion + species1_proportion >= 0.85) && (species0_proportion >= 0.15) && (species1_proportion >= 0.15)
                        outcome = "two_species"
                    else
                        outcome = "blended"
                    end
                    HI_NL_all_inds = [calc_traits_additive(genotypes_F[:, neutral_loci, :]); calc_traits_additive(genotypes_M[:, neutral_loci, :])]
                    species0_proportion_NL = sum(HI_NL_all_inds .== 0) / length(HI_NL_all_inds)
                    species1_proportion_NL = sum(HI_NL_all_inds .== 1) / length(HI_NL_all_inds)
                    if save_each_sim
                        @save string("HZAM_Sym_Julia_results_GitIgnore/simulation_data.", run_name, ".jld2") outcome functional_HI_all_inds species0_proportion species1_proportion HI_NL_all_inds species0_proportion_NL species1_proportion_NL
                    end
                end
                println(run_name, "  outcome was: ", outcome)
                outcome_array[i, j, k] = outcome
            end # of S_AM loop
        end # of w_hyb loop   
    end # of replicate loop

    if save_outcomes_JL
        filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set",set_name,"_surv",short_survFitnessMethod,"_ecolDiff",ecolDiff,"_growthrate",intrinsic_R,"_K",K_total,"_FL",total_functional_loci,"_NL",num_neutral_loci,"_gen",max_generations,"_SC",per_reject_cost,".jld2")
        save_object(filename, outcome_array)
    end
 
    if save_outcomes_csv
        for i in 1:size(outcome_array, 3)
            filename = string("HZAM_Sym_Julia_results_GitIgnore/outcomeArray_set",set_name,"_surv",short_survFitnessMethod,"_ecolDiff",ecolDiff,"_growthrate",intrinsic_R,"_K",K_total,"_FL",total_functional_loci,"_NL",num_neutral_loci,"_gen",max_generations,"_SC",per_reject_cost,"_rep",replications[i])
            CSV.write(filename, Tables.table(outcome_array[:,:,i]), writeheader=false)
        end 
    end
    return outcome_array
end



#### functions for summarizing and plotting results

# convert outcome array (an array of strings) to categorical array
function convert_to_cat_array(outcome_array)
    cat_outcome_array = compress(CategoricalArray(outcome_array))
    levels!(cat_outcome_array, ["extinction", "blended", "one_species", "two_species"])
    return cat_outcome_array
end 

# function for plotting grid of pie charts showing distribution of four outcomes, using categorical array as input
function plot_all_outcomes(cat_outcome_array)
    num_outcome_types = length(levels(cat_outcome_array))
    outcome_counts = Array{Int16, 2}(undef, num_outcome_types, (size(cat_outcome_array, 1)*size(cat_outcome_array, 2) )) 
    for i in 1:size(cat_outcome_array, 1) 
        for j in 1:size(cat_outcome_array, 2) 
            for outcome_num in 1:num_outcome_types
                outcome_counts[outcome_num, j + (i-1)*size(cat_outcome_array, 2)] = sum(cat_outcome_array[i,j,:] .== levels(cat_outcome_array)[outcome_num])
            end
        end
    end
    colors_of_outcomes = [RGB(0,0,0), plasma[0.525], plasma[0.2], plasma[0.9]] # colors for 4 outcome categories
    pie(outcome_counts, layout = grid(size(cat_outcome_array, 1), size(cat_outcome_array, 2)), legend = false, palette = colors_of_outcomes, margin = -2.0mm)
    plot!(size=(800,1300))
end

function get_most_common_outcomes(cat_outcome_array)
    levels_of_outcomes = levels(cat_outcome_array)
    num_outcome_types = length(levels_of_outcomes)
    most_common_outcomes = CategoricalArray{String, 2}(undef, size(cat_outcome_array, 1), size(cat_outcome_array, 2))
    levels!(most_common_outcomes, ["extinction", "blended", "one_species", "two_species"])
    for i in 1:size(cat_outcome_array, 1) 
        for j in 1:size(cat_outcome_array, 2)
            outcome_counts = Vector{Int}(undef, num_outcome_types) 
            for outcome_num in 1:num_outcome_types
                outcome_counts[outcome_num] = sum(cat_outcome_array[i,j,:] .== levels(cat_outcome_array)[outcome_num])
            end
            outcomes_with_max_count = findall(outcome_counts .== maximum(outcome_counts))
            if length(outcomes_with_max_count) == 1
                most_common_outcomes[i,j] = levels_of_outcomes[outcomes_with_max_count[1]]
            elseif length(outcomes_with_max_count) >= 2
                # if tie in outcome count, choose one randomly
                most_common_outcomes[i,j] = levels_of_outcomes[sample(outcomes_with_max_count)]
            end
        end
    end
    return most_common_outcomes
end

function plot_common_outcomes(common_outcome_array)
    # make heat map of outcomes
    colors_of_outcomes = [RGB(0,0,0), plasma[0.525], plasma[0.2], plasma[0.9]] # colors for 4 outcome categories
    one_outcome_array = reverse(common_outcome_array, dims=1)
    x_midpoints = log10.([1, 3, 10, 30, 100, 300, 1000, 5000])  # the S_AM values, with Inf convert to 3000 for graphing 
    w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # for just one run, just put one number in this and next line
    y_midpoints = reverse(w_hyb_set)
    min_color = minimum(levelcode.(one_outcome_array)) # this and next line needed to choose the proper colors for the figure
    max_color = maximum(levelcode.(one_outcome_array))
    heatmap(x_midpoints, y_midpoints, one_outcome_array, c = colors_of_outcomes[min_color:max_color], yflip = false, tick_direction = :out, colorbar = false, size = (440,310), framestyle = :box)
    xaxis!("Strength of conspecific mate preference")
    xticklabels = ["1", "3", "10", "30", "100", "300", "1000", "complete"]
    plot!(xticks=(x_midpoints, xticklabels))
    yaxis!("Hybrid fitness")
    yticklabels = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "", "", "1.0"]
    plot!(yticks=(y_midpoints, yticklabels), tick_direction = :out)
    # add white lines
    plot!([xlims()[1], xlims()[2]], [0.05, 0.05], linecolor = :white, widen = false, legend = false, linewidth=3)
    x_for_line = mean(x_midpoints[[length(x_midpoints)-1 length(x_midpoints)]])
    plot!([x_for_line, x_for_line], [ylims()[1], ylims()[2]], linecolor = :white, widen = false, legend = false, linewidth=3)
end

function make_and_save_figs(ResultsFolder, RunName, RunOutcomes)
    cat_RunOutcomes = convert_to_cat_array(RunOutcomes)
    display(plot_all_outcomes(cat_RunOutcomes))
    savefig(string(ResultsFolder,"/",RunName,"_AllOutcomes.png"))
    savefig(string(ResultsFolder,"/",RunName,"_AllOutcomes.pdf"))
    most_common_outcomes = get_most_common_outcomes(cat_RunOutcomes)
    display(plot_common_outcomes(most_common_outcomes))
    savefig(string(ResultsFolder,"/",RunName,"_MostCommonOutcomes.png"))
    savefig(string(ResultsFolder,"/",RunName,"_MostCommonOutcomes.pdf"))
end


