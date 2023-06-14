include("src/simulation.jl")

run_one_HZAM_sim(0.7, 10, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=10000, max_generations=1000,
    sigma_disp=0.05, sympatry=false,
    sigma_comp=0.1, do_plot=true, plot_int=1, optimize=false)