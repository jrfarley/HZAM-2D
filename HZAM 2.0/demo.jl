include("src/simulation.jl")

K = Integer(trunc(1000*(1/sqrt(2*pi*0.01^2)))) # Number of individuals so that the density is the same as in 1D model with 1000 individuals

K = 40342

println(K)

run_one_HZAM_sim(1.0, 10, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=50,
    sigma_disp=0.02, sigma_comp=0.01, do_plot=true, plot_int=1)