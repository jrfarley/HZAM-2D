include("HZAM/src/HZAM.jl")

import .HZAM

K = 20000

println(K)

HZAM.run_one_HZAM_sim(0.8, 1, 1, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=1000,
    sigma_disp=0.07, sigma_comp=0.01, do_plot=true, plot_int=1, total_loci=10)
