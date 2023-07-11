@testset "initialize_population_1" begin
    n = Population.NUM_DEMES
    K_total = 2 * (n^2)
    ecolDiff = 0
    total_loci = 3
    intrinsic_R = 1.1
    sigma_comp = 0.01

    pd = PopulationData(K_total,
        ecolDiff,
        total_loci,
        intrinsic_R,
        sigma_comp)

        genotypes_F = vcat([d.genotypes_F for d in pd.population]...)
        genotypes_M = vcat([d.genotypes_M for d in pd.population]...)
        locations_F = vcat([d.locations_F for d in pd.population]...)
        locations_M = vcat([d.locations_M for d in pd.population]...)
            
    

    @test count(i->i==[0 0 0; 0 0 0], genotypes_F) == n^2 / 2
    @test count(i->i==[0 0 0; 0 0 0], genotypes_M) == n^2 / 2
    @test count(i->i==[1 1 1; 1 1 1], genotypes_F) == n ^2 / 2
    @test count(i->i==[1 1 1; 1 1 1], genotypes_M) == n^2 /2 
    @test length(genotypes_F) == n^2
    @test length(genotypes_M) == n^2 
    @test length(locations_F) == n^2
    @test length(locations_M) == n^2

    [@test length(d.locations_F) == 1 for d in pd.population]
    [@test length(d.locations_M) == 1 for d in pd.population]
    [@test length(d.genotypes_F) == 1 for d in pd.population]
    [@test length(d.genotypes_M) == 1 for d in pd.population]
    [@test assign_zone(pd.population[i].locations_F[1])==i for i in eachindex(IndexCartesian(), pd.population)]
    [@test assign_zone(pd.population[i].locations_M[1])==i for i in eachindex(IndexCartesian(), pd.population)]
    [@test pd.population[i].mitochondria_F[1]==0 && i[1]<=5 ||pd.population[i].mitochondria_F[1]==1 && i[1]>5 for i in eachindex(IndexCartesian(), pd.population)]
    [@test pd.population[i].mitochondria_M[1]==0 && i[1]<=5 ||pd.population[i].mitochondria_M[1]==1 && i[1]>5 for i in eachindex(IndexCartesian(), pd.population)]
    [@test pd.population[i].genotypes_F[1]==[0 0 0; 0 0 0] && i[1]<=5 ||pd.population[i].genotypes_F[1]==[1 1 1; 1 1 1] && i[1]>5 for i in eachindex(IndexCartesian(), pd.population)]
    [@test pd.population[i].genotypes_M[1]==[0 0 0; 0 0 0] && i[1]<=5 ||pd.population[i].genotypes_M[1]==[1 1 1; 1 1 1] && i[1]>5 for i in eachindex(IndexCartesian(), pd.population)]

    [@test length(pd.growth_rates_F[i])==1 for i in eachindex(IndexCartesian(), pd.population)]
end

@testset "calc_traits_additive" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

    @test calc_traits_additive(genotypes, 1:3) == [0, 1]

    genotypes = [[1 0 0; 0 0 0], [0 1 0; 1 0 0], [1 1 1; 1 1 1]]

    @test calc_traits_additive(genotypes, 1:2) == [0.25, 0.5, 1]
end

@testset "disperse_individual" begin
    locations = Location[]
    geographic_limits = [Location(0.0f0, 0.0f0), Location(1.0f0, 1.0f0)]
    sigma_disp = 0.01

    num_fail = 0

    for i in 1:1000
        new_location = Location(Location(0.5f0, 0.5f0), sigma_disp, geographic_limits)
        push!(locations, new_location)
        if new_location.x < 0.5 - 2 * sigma_disp || new_location.x > 0.5 + 2 * sigma_disp || new_location.y < 0.5 - 2 * sigma_disp || new_location.y > 0.5 + 2 * sigma_disp
            num_fail += 1
        end
    end

    average_x = sum([l.x for l in locations]) / 1000

    average_y = sum([l.y for l in locations]) / 1000

    @test num_fail < 70

    @test average_x > 0.495 && average_x < 0.505
    @test average_y > 0.495 && average_y < 0.505

end

@testset "mean" begin
    genotypes = [[1 0 1; 0 1 0], [0 0 1; 1 0 0], [1 1 1; 1 1 1]]

    @test Population.mean(genotypes[1], collect(1:2)) == 0.5

    @test Population.mean(genotypes[2], [2]) == 0
    @test Population.mean(genotypes[3], collect(1:3)) == 1
end

@testset "assign_zone" begin
    locations = [Location(0.01f0, 0f0), Location(0.32f0, 0.75f0), Location(0.5f0, 0.5f0)]

    zones = assign_zone.(locations)

    @test zones == [CartesianIndex(1,1), CartesianIndex(4,8), CartesianIndex(6,6)]
end