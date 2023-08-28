# Population
Types and functions used for managing population data (locations, genotypes, and growth rates).

## Types
```@docs
HZAM.Population.Location
HZAM.Population.Zone
HZAM.Population.PopulationData
```

## Functions
```@docs
HZAM.Population.assign_zone(location::HZAM.Population.Location)
HZAM.Population.calc_ind_useResourceA(competition_traits::Vector, ecolDiff::Real)
HZAM.Population.calc_ind_useResourceB(competition_traits::Vector, ecolDiff::Real)
HZAM.Population.max_radius_squared(location::HZAM.Population.Location, t::Real, max_dist::Real)
HZAM.Population.calc_ideal_densities(
    K_total::Integer,
    sigma_comp::Real,
    locations_F::Vector{<:HZAM.Population.Location},
    max_dist::Real
)
HZAM.Population.calc_squared_distances(location_list::Vector{<:HZAM.Population.Location}, focal_location::HZAM.Population.Location)
HZAM.Population.get_surrounding_zones(zone_index::CartesianIndex)
HZAM.Population.calc_real_densities(
    neighbourhood::Matrix{HZAM.Population.Zone},
    locations_F::Vector{<:HZAM.Population.Location},
    max_dist::Real,
    sigma_comp::Real
)
HZAM.Population.calc_growth_rates(
    population::Matrix{HZAM.Population.Zone},
    zone_index::CartesianIndex,
    K_total::Integer,
    sigma_comp::Real,
    intrinsic_R::Real
)
HZAM.Population.genotype_mean(
    genotype::Matrix{<:Integer},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
HZAM.Population.calc_survival_fitness_hetdisadvantage(
    genotype::Matrix{<:Integer},
    w_hyb::Real
)
HZAM.Population.calc_survival_fitness_epistasis(genotype::Matrix{<:Integer},
    hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
    w_hyb::Real,
    beta::Real=1
)
```