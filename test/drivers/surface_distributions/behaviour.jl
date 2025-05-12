@testset "behaviour - rand, cdf, pdf" begin

    # rand - correct surface radius
    @test all(r -> rand(EqualSurfaceDistribution(r))[1] == r, rand(1000))
    @test all(r -> rand(SolarSurfaceDistribution(r))[1] == r, rand(1000))

    # pdf - equal distribution is equal at equal latitudes
    lat = rand()*pi-pi/2
    @test pdf(EqualSurfaceDistribution(100), (1, rand()*2pi-pi, lat)) == 
          pdf(EqualSurfaceDistribution(100), (1, rand()*2pi-pi, lat))

    # pdf - solar surface distribution zero on night-side
    @test pdf(SolarSurfaceDistribution(100), (1, -pi/2-rand(), rand()*pi-pi/2)) == 0
    @test pdf(SolarSurfaceDistribution(100), (1,  pi/2+rand(), rand()*pi-pi/2)) == 0

    # pdf behavior (higher latitude, lower likelihood)
    @test all(lat -> pdf(EqualSurfaceDistribution(1), (1, 0, lat)) > pdf(EqualSurfaceDistribution(1), (1, 0, lat+1e-3)), rand(1000)*(pi/2 - 1e-3))
    @test all(lat -> pdf(SolarSurfaceDistribution(1), (1, 0, lat)) > pdf(SolarSurfaceDistribution(1), (1, 0, lat+1e-3)), rand(1000)*(pi/2 - 1e-3))

    # pdf behavior - SolarSurfaceDistribution - higher absolute longitude, lower likelihood
    @test all(lon -> pdf(SolarSurfaceDistribution(1), (1, lon, 0)) > pdf(SolarSurfaceDistribution(1), (1, lon+1e-3, 0)), rand(1000)*pi/2)

    # cdf behavior (full domain)
    @test isapprox(cdf(EqualSurfaceDistribution(100), (1, -pi, -pi/2), (1, pi, pi/2)), 1; rtol=1e-4)
    @test isapprox(cdf(SolarSurfaceDistribution(100), (1, -pi/2, -pi/2), (1, pi/2, pi/2)), 1; rtol=1e-4)

    # cdf behavior (empty domain)
    x = rand(3)
    @test isapprox(cdf(EqualSurfaceDistribution(100), x, x), 0; rtol=1e-4)
    @test isapprox(cdf(SolarSurfaceDistribution(100), x, x), 0; rtol=1e-4)
end