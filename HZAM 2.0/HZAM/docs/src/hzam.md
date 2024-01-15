# HZAM

## Index
```@index
Pages = ["population.md", "mating.md", "data_analysis.md", "plot_data.md", "index.md"]
```

## Functions
```@docs
HZAM.run_one_HZAM_sim(
        w_hyb::Real, 
        S_AM::Real,
        intrinsic_R::Real;
    )
HZAM.run_HZAM_set(set_name::String, intrinsic_R::Real, ecolDiff::Real;)
HZAM.load_from_folder(dir::String)
HZAM.load_from_csv(filepath::String)
HZAM.convert_to_CSV(
    outcome_array::Array{<:Real},
    w_hyb_array::Array{<:Real},
    S_AM_array::Array{<:Real},
    name::String
)
HZAM.plot_fitnesses(fitnesses::Vector{<:Dict})
HZAM.plot_population_tracking_data(filepath::String)
HZAM.summarize_gene_correlations(dir::String)
HZAM.plot_phenotypes(phenotypes;)
HZAM.set_results_folder(dir::String)
```
