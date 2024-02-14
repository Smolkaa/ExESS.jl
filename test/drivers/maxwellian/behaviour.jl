@testset "behaviour - rand, cdf, pdf" begin

    # boundary behavior (zero temperature => zero speed/velocity)
    @test rand(MBSpeedDistribution(0, amu2kg(rand()*100)))        == 0
    @test rand(MBFluxSpeedDistribution(0, amu2kg(rand()*100)))    == 0
    @test rand(MBVelocityDistribution(0, amu2kg(rand()*100)))     == (0,0,0)
    @test rand(MBFluxVelocityDistribution(0, amu2kg(rand()*100))) == (0,0,0)

    # boundary behavior (zero temperature => zero probability)
    @test pdf(MBSpeedDistribution(0, amu2kg(rand()*100)), rand()*1000)        == 0
    @test pdf(MBFluxSpeedDistribution(0, amu2kg(rand()*100)), rand()*1000)    == 0

    # boundary behavior (zero mass => infinite speed)
    @test isinf(rand(MBSpeedDistribution(rand()*400, 0)))
    @test isinf(rand(MBFluxSpeedDistribution(rand()*400, 0)))

    # boundary behavior (zero mass => zero probability)
    @test pdf(MBSpeedDistribution(rand()*400, 0), rand()*1000)     == 0
    @test pdf(MBFluxSpeedDistribution(rand()*400, 0), rand()*1000) == 0

    #::. cdf behavior (full domain)
    @test isapprox(cdf(MBAzimuthDistribution(rand()*400, amu2kg(rand()*100)), -pi, pi), 1; rtol=1e-4)
    @test isapprox(cdf(MBFluxAzimuthDistribution(rand()*400, amu2kg(rand()*100)), -pi, pi), 1; rtol=1e-4)
    @test isapprox(cdf(MBElevationDistribution(rand()*400, amu2kg(rand()*100)), -pi/2, pi/2), 1; rtol=1e-4)
    @test isapprox(cdf(MBFluxElevationDistribution(rand()*400, amu2kg(rand()*100)), 0, pi/2), 1; rtol=1e-4) # half domain for MBF
    @test isapprox(cdf(MBSpeedDistribution(rand()*400, amu2kg(rand()*100)), 0, 100_000), 1; atol=1e-4)
    @test isapprox(cdf(MBFluxSpeedDistribution(rand()*400, amu2kg(rand()*100)), 0, 100_000), 1; atol=1e-4)
    
    #::. cdf behavior (empty domain)
    @test isapprox(cdf(MBAzimuthDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; rtol=1e-4)
    @test isapprox(cdf(MBFluxAzimuthDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; rtol=1e-4)
    @test isapprox(cdf(MBElevationDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; rtol=1e-4)
    @test isapprox(cdf(MBFluxElevationDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; rtol=1e-4)
    @test isapprox(cdf(MBSpeedDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; atol=1e-4)
    @test isapprox(cdf(MBFluxSpeedDistribution(rand()*400, amu2kg(rand()*100)), rand().*(1,1) ...), 0; atol=1e-4)

end