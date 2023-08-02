include("HZAM/src/HZAM.jl")
import .HZAM
import HypothesisTests.EqualVarianceTTest

K = 20000

loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15, 16], [17, 18, 19, 20], 1:4, 5:8, 9:12, 13:16))
#=
outcome, pd = HZAM.run_one_HZAM_sim(0.8, 100, 1, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=1000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true,      plot_int=10,
    total_loci=20,
    female_mating_trait_loci=1:4,
    male_mating_trait_loci=5:8,
    competition_trait_loci=9:12,
    hybrid_survival_loci=13:16)


=#
filepath = "HZAM_Sym_Julia_results_GitIgnore/simulation_outcomes/genotypes_ecolDiff0_w_hyb0.4_S_AM1.jld2"

genotypes = HZAM.load_genotypes(filepath)

function calc_chi_squared(genotypes, locus)
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

function t_test(genotypes, loci)
    function calc_heterozygosity(locus)
        count(x->x==[0; 1] || x == [1; 0], [g[:, locus] for g in genotypes]) / length(genotypes)
    end


    heterozygosities = calc_heterozygosity.(loci.overall)

    neutral_heterozygosities = heterozygosities[loci.neutral]

    for trait in keys(loci)[3:7]
        println(trait)
        println(EqualVarianceTTest(heterozygosities[loci.trait], neutral_heterozygosities))
    end
end

function t_test_chi_squared(genotypes, loci)
    chi_squareds = calc_chi_squared.(Ref(genotypes), loci.overall)

    for trait in keys(loci)[3:7]
        println(trait)
        println(EqualVarianceTTest(chi_squareds[getfield(loci,trait)], chi_squareds[loci.neutral]))
    end
end

t_test_chi_squared(genotypes, loci)
HZAM.chi_squared_traits(genotypes, loci)
println(HZAM.find_extinct_alleles(genotypes))
println(calc_chi_squared.(Ref(genotypes), loci.neutral))
HZAM.create_gene_plot(genotypes, loci, 1000)