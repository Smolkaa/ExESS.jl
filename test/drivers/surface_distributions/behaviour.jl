@testset "behaviour - rand, cdf, pdf" begin

    # correct surface radius
    r = rand()*100
    @test rand(EqualSurfaceDistribution(r))[1] == r
    @test rand(SolarSurfaceDistribution(r))[1] == r

    # equal distribution is equal at equal latitudes
    lat = rand()*pi-pi/2
    @test pdf(EqualSurfaceDistribution(100), (1, rand()*2pi-pi, lat)) == 
          pdf(EqualSurfaceDistribution(100), (1, rand()*2pi-pi, lat))

    # solar surface distribution zero on night-side
    @test pdf(SolarSurfaceDistribution(100), (1, -pi/2-rand(), rand()*pi-pi/2)) == 0
    @test pdf(SolarSurfaceDistribution(100), (1,  pi/2+rand(), rand()*pi-pi/2)) == 0

    # cdf behavior (full domain)
    @test isapprox(cdf(EqualSurfaceDistribution(100), (1, -pi, -pi/2), (1, pi, pi/2)), 1; rtol=1e-4)
    @test isapprox(cdf(SolarSurfaceDistribution(100), (1, -pi/2, -pi/2), (1, pi/2, pi/2)), 1; rtol=1e-4)

    # cdf behavior (empty domain)
    x = rand(3)
    @test isapprox(cdf(EqualSurfaceDistribution(100), x, x), 0; rtol=1e-4)
    @test isapprox(cdf(SolarSurfaceDistribution(100), x, x), 0; rtol=1e-4)
end