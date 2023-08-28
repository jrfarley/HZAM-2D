# PlotData
Functions for plotting data while the simulation is running.

## Functions
```@docs
HZAM.PlotData.create_population_plot(
    hybrid_indices_functional::Vector,
    locations::Vector,
    save_plot::Bool
)
HZAM.PlotData.update_population_plot(
    hybrid_indices_functional::Vector,
    locations::Vector,
    generation::Integer,
    save_plot::Bool
)
HZAM.PlotData.scale_curve(curve, Index::Integer)
HZAM.PlotData.count_genotypes(
        range::Union{UnitRange{<:Integer},Vector{<:Integer}},
        genotypes::Vector{<:Matrix{<:Integer}}
    )
HZAM.PlotData.get_expected(x::Real, n::Integer)
HZAM.PlotData.create_gene_plot(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple,
    save_plot::Bool
)
HZAM.PlotData.update_gene_plot(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::NamedTuple,
        generation::Integer,
        save_plot::Bool
    )
```