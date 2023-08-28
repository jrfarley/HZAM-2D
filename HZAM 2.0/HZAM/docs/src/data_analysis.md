# DataAnalysis
Types and functions used for analyzing data produced by the simulation.

## Types
```@docs
HZAM.DataAnalysis.SimParams
HZAM.DataAnalysis.PopulationTrackingData
HZAM.DataAnalysis.SpatialData
HZAM.DataAnalysis.OutputData
```

## Functions

```@docs
HZAM.DataAnalysis.calc_position(sigmoid_curves::Vector{<:Vector{<:Real}})
HZAM.DataAnalysis.calc_variance(positions::Vector{<:Real})
HZAM.DataAnalysis.calc_sigmoid_curve(locations_x::Vector{<:Real}, hybrid_indices::Vector{<:Real})
HZAM.DataAnalysis.calc_sigmoid_curves(locations::Vector, hybrid_indices::Vector{<:Real})
HZAM.DataAnalysis.sigmoid(x::Vector{<:Real}, p::Vector{<:Real})
HZAM.DataAnalysis.calc_width(sigmoid_curve::Vector{<:Real})
HZAM.DataAnalysis.average_width(sigmoid_curves::Vector{<:Vector{<:Real}})
HZAM.DataAnalysis.calc_length(sigmoid_curves::Vector{<:Vector{<:Real}})
HZAM.DataAnalysis.calc_all_gene_flow(
    genotypes::Vector{<:Matrix{<:Real}},
    hybrid_indices_functional::Vector{<:Real},
    loci::NamedTuple
)
HZAM.DataAnalysis.calc_all_cline_widths_and_positions(
        genotypes::Vector{<:Matrix{<:Real}},
        locations::Vector,
        loci::NamedTuple
    )
HZAM.DataAnalysis.calc_overlap_in_range(
        locations_x::Vector{<:Real},
        hybrid_indices_functional::Vector{<:Real}
    )
HZAM.DataAnalysis.calc_overlap_overall(
        locations_x::Vector,
        hybrid_indices_functional::Vector{<:Real},
        sorted_indices::Vector{<:Vector{<:Integer}}
    )
HZAM.DataAnalysis.calc_bimodality_in_range(
    sigmoid_curve::Vector{<:Real},
    locations_x::Vector{<:Real},
    hybrid_indices::Vector{<:Real},
    sigma_disp::Real
)
HZAM.DataAnalysis.calc_bimodality_overall(
        sigmoid_curves::Vector{<:Vector{<:Real}},
        sorted_indices::Vector{<:Vector{<:Integer}},
        locations_x::Vector{<:Real},
        hybrid_indices::Vector{<:Real},
        sigma_disp::Real
    )
HZAM.DataAnalysis.calc_linkage_diseq(
        genotypes::Vector{<:Matrix{<:Integer}},
        l1::Integer,
        l2::Integer
    )
HZAM.DataAnalysis.calc_average_linkage_diseq(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple
)
HZAM.DataAnalysis.calc_trait_correlation(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci_range1::Union{UnitRange{<:Integer},Vector{<:Integer}},
        loci_range2::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )
HZAM.DataAnalysis.calc_all_trait_correlations(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple
)
HZAM.DataAnalysis.sort_y(y_locations::Vector{<:Real})
HZAM.DataAnalysis.sort_locations(A::AbstractArray{<:Real}, bin_size::Real)
HZAM.DataAnalysis.average_gene_data(gene_data::Vector{<:NamedTuple})
HZAM.DataAnalysis.calc_traits_additive(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
HZAM.DataAnalysis.average_data_per_phenotype(
        data::Vector{<:Integer},
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )
HZAM.DataAnalysis.count_phenotypes_at_loci(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
HZAM.DataAnalysis.find_fixed_alleles(genotypes::Vector{<:Matrix{<:Integer}})
HZAM.DataAnalysis.calc_chi_squared(genotypes::Vector{<:Matrix{<:Integer}}, locus::Integer)
```