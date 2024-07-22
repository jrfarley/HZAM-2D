"Types and functions used for analyzing data produced by the simulation."
module DataAnalysis

using LsqFit: curve_fit

mean(itr) = sum(itr) / length(itr)

include("used_data_analysis.jl")
include("unused_data_analysis.jl")
end