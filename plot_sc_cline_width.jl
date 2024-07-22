include("HZAM 2.0/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using Plots
using Statistics

outcome_array = HZAM.load_from_folder2("HZAM-J_2D_results/search_cost_cline_sc_0.05")
widths = [A.hybrid_zone_width[21:40] for (A, S_AM) in outcome_array]

mid_sc = [mean(w) for w in widths]
S_AMs = [S_AM for (A, S_AM) in outcome_array]
σs = [std(w) for w in widths]

perm = sortperm(S_AMs)

plot(S_AMs[perm], mid_sc[perm], yerror=σs[perm], xaxis = :log, label="Search cost: 0.05")

outcome_array = HZAM.load_from_folder2("HZAM-J_2D_results/search_cost_cline_sc_0.01")
widths = [A.hybrid_zone_width[21:40] for (A, S_AM) in outcome_array]

weak_sc = [mean(w) for w in widths]
S_AMs = [S_AM for (A, S_AM) in outcome_array]
σs = [std(w) for w in widths]

perm = sortperm(S_AMs)

plot!(S_AMs[perm], weak_sc[perm], yerror=σs[perm], xaxis = :log, label="Search cost: 0.01")

outcome_array = HZAM.load_from_folder2("HZAM-J_2D_results/search_cost_cline_sc_0.1")
widths = [A.hybrid_zone_width[21:40] for (A, S_AM) in outcome_array]

strong_sc = [mean(w) for w in widths]
S_AMs = [S_AM for (A, S_AM) in outcome_array]
σs = [std(w) for w in widths]

perm = sortperm(S_AMs)

plot!(S_AMs[perm], strong_sc[perm], yerror=σs[perm], xaxis = :log, label="Search cost: 0.1")

gui()
readline()
