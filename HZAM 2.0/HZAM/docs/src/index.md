# HZAM

## Contents
```@contents
```

## Functions
```@docs
HZAM.run_one_HZAM_sim(
        w_hyb::Real, 
        S_AM::Real, 
        ecolDiff::Real, 
        intrinsic_R::Real;
    )
HZAM.run_HZAM_set(set_name::String, intrinsic_R::Real, ecolDiff::Real;)
HZAM.plot_output_field(
    outcomes::Array{:Real},
    sim_params::Array{<:HZAM.DataAnalysis.SimParams}
)
HZAM.load_from_folder(dir::String)
HZAM.load_from_csv(filepath::String)
HZAM.convert_to_CSV(
    outcome_array::Array{<:HZAM.DataAnalysis.OutputData},
    field_name::Symbol,
    output_folder::String
)
HZAM.plot_fitnesses(fitnesses::Vector{<:Dict})
HZAM.plot_population_tracking_data(filepath::String)
HZAM.summarize_gene_correlations(dir::String)
HZAM.plot_phenotypes(phenotypes;)
```

## Index
```@index
```