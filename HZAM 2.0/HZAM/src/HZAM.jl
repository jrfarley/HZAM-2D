#=import Pkg
Pkg.activate("./HZAM")
=#

module HZAM
include("data_analysis.jl") # functions for analyzing data (calculating width/length of hybrid zone, gene flow, bimodality, etc.)
include("population.jl") # data types and functions for initializing population with genotypes, locations, growth rates, etc.
include("mating.jl") # functions for finding a mate and determining match strength
include("plot_data.jl")

using .DataAnalysis
using .Population
using .Mating
using .Plot_Data

export DataAnalysis, Population, Mating, demo # export modules for testing

global demo = 10

include("simulation.jl") # main function for running a simulation

end # module HZAM
