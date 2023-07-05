include("src/simulation.jl")

#K = Integer(trunc(1000*(1/sqrt(2*pi*0.01^2)))) # Number of individuals so that the density is the same as in 1D model with 1000 individuals

# K should be 40,002 individuals to give the same density as the 1D model

# If there are fewer than a few thousand individuals the population will go extinct

K = 5000

println(K)

run_one_HZAM_sim(0.9, 1, 0, 1.1; # these values are 
    # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
    K_total=K, max_generations=10000,
    sigma_disp=0.05, sigma_comp=0.01, do_plot=true, plot_int=1)