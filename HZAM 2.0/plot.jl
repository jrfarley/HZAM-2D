
using GLMakie
using QuadGK
# initialize plot

fig = Figure(resolution=(3840, 2160))

# add axis

ax1 = fig[1, 1] = Axis(fig,
    title="Fitness by phenotype",
    titlegap=48, titlesize=60,
    xautolimitmargin=(0, 0), xgridwidth=2, xticklabelsize=36,
    xticks=LinearTicks(20), xticksize=18,
    yautolimitmargin=(0, 0), ygridwidth=2, yticklabelpad=14,
    yticklabelsize=36, yticks=LinearTicks(20), yticksize=18, limits=(-3, 3, -0.1, 1.1),
    xlabel="X Location", ylabel="Fitness"
)

vlines!(ax1, [0], linewidth=2)
hlines!(ax1, [0], linewidth=2)


sg = SliderGrid(
    fig[1, 2],
    (label="S_AM", range=2:10:500, format="{:.1f}", startvalue=50),
    (label="w_hyb", range=0:0.01:1, format="{:.1f}", startvalue=0.7),
    (label="σ", range=0.1:0.1:1, format="{:.1f}", startvalue=0.1),
    (label="mating_trait", range=0:0.1:1, format="{:.1f}", startvalue=0.75),
    width=350,
    tellheight=false
)

S_AM = sg.sliders[1].value # assortative mating strength

w_hyb = sg.sliders[2].value # hybrid fitness

σ = sg.sliders[3].value # standard deviation for the normal distribution governing the probability of two individuals meeting

trait = sg.sliders[4].value # hybrid index for the green line on the plot

# width of female acceptance curve for male trait
function pref_SD(S_AM)::Float32
    sqrt(-(1 / (2 * log(1 / S_AM))))
end

# the average hybrid index at any point along the x axis
function gene_distribution(x)
    1 / (1 + exp(-10 * x))
end

#  generic normal distribution function
function normal(x, μ, σ)
    (1 / (σ * sqrt(2 * π))) * exp(-(1 / 2) * ((x - μ) / σ)^2)
end

# probability of female accepting mate
function mating_prob(x, μ, pref_SD)
    exp(-((x - μ)^2) / (2 * pref_SD))
end

# probability of offspring surviving
function survival(trait, w_hyb)
    1 - (1 - w_hyb) * 4 * trait * (1 - trait)
end

# calculate the fitness of a male at location x with a given hybrid index with respect to the female at location μ
function fitness(trait, x, μ, σ, pref_SD, w_hyb)
    mating_prob(trait, gene_distribution(μ), pref_SD) *
    normal(x, μ, σ) * survival((trait + gene_distribution(μ)) / 2, w_hyb)
end

# calculate the overall fitness of a male by integrating the fitness function over the whole range
function integrate_fitness(trait, x, σ, pref_SD, w_hyb)
    quadgk(t -> fitness(trait, x, t, σ, pref_SD, w_hyb), -15, 15)[1]
end

x = -10:0.01:10
trait0 = @lift(integrate_fitness.(Ref(0), x, Ref($σ), Ref(pref_SD.(Ref($S_AM))), Ref($w_hyb)))
trait1 = @lift(integrate_fitness.(Ref(1), x, Ref($σ), Ref(pref_SD.(Ref($S_AM))), Ref($w_hyb)))
trait2 = @lift(integrate_fitness.(Ref(0.5), x, Ref($σ), Ref(pref_SD.(Ref($S_AM))), Ref($w_hyb)))
trait3 = @lift(integrate_fitness.(Ref($trait), x, Ref($σ), Ref(pref_SD.(Ref($S_AM))), Ref($w_hyb)))



line1 = lines!(ax1, x, trait0, color=:blue, linewidth=5)
line2 = lines!(ax1, x, trait1, color=:red, linewidth=5)
line3 = lines!(ax1, x, trait2, color=:purple, linewidth=5)
line4 = lines!(ax1, x, trait3, color=:green, linewidth=5)

axislegend(ax1,
    [line1, line2, line3, line4],
    ["HI=0", "HI=1", "HI=0.5", "HI=mating_trait"],
    "Legend", position=:rt)



display(fig)
readline()