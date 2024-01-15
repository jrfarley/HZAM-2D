"Types and functions used for analyzing data produced by the simulation."
module DataAnalysis
export OutputData, SimParams, PopulationTrackingData

using LsqFit: curve_fit
using Statistics: mean

import HypothesisTests.EqualVarianceTTest

include("used_data_analysis.jl")
#include("unused_data_analysis.jl")
end