using Colors, ColorSchemes
import ColorSchemes.plasma
using GLMakie
using JLD2

filenames = [
    "w_hyb0_S_AM1_ecolDiff1_full_pleiotropy", # narrow tension zone
    "w_hyb0.8_S_AM1_ecolDiff1_full_pleiotropy", # narrow tension zone
    "w_hyb0.95_S_AM1_no_pleiotropy", # broad tension zone
    "w_hyb0.8_S_AM100_full_pleiotropy", # hybridization and overlap
    "w_hyb0_S_AM100_no_pleiotropy", #limited overlap 
    "w_hyb1_S_AM_inf_ecolDiff1_full_pleiotropy", #extensive overlap
    "w_hyb1_S_AM1_ecolDiff1_full_pleiotropy", # blended 3
    "w_hyb1_S_AM100_ecolDiff1_magic_cue" # blended 2
]

plotnames = ["Narrow Tension Zone", "Narrow Tension Zone", "Broad Tension Zone", "Hybridization and Overlap", "Limited Overlap", "Extensive Overlap", "Blended", "Blended", "Blended"]

function calculate_output(hybrid_indices, locations)
    output = Matrix{Real}(undef, 7, 50)
    x_bins = collect(0.02:0.02:1)
    println(length(x_bins))
    for i in 1:50
        for j in 0:6
            indices = findall(x -> x_bins[i] - 0.02 <= x < x_bins[i], locations)
            if length(indices) ≠ 0
                output[j+1, i] = count(ind -> hybrid_indices[ind] ≈ j / 6, indices) / length(indices)
            else
                output[j+1, i] = 0
            end
        end
    end

    return output
end

function make_plot()
    f = Figure()
    ax = []
    for i in 1:3
        for j in 1:3
            push!(ax, Axis(f[i, j], xlabel="generation", ylabel="phenotype", title=plotnames[(i-1)*3+j]))
        end
    end

    for i in 1:length(filenames)
        @load filenames[i] hybrid_indices locations
        output = calculate_output(hybrid_indices, locations)
        heatmap!(ax[i], 0:0.02:1, 0:6, output', colormap=:Greys)
    end
    display(f)
    readline()
end


make_plot()