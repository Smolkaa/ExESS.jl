print("TESTING: mercury > constants.jl")

@testset verbose=true "constants.jl ............................." begin
    #::. existence
    @test isdefined(ExESS, :MERCURY_DAY)
    @test isdefined(ExESS, :MERCURY_MASS)
    @test isdefined(ExESS, :MERCURY_ORBITAL_ECCENTRICITY)
    @test isdefined(ExESS, :MERCURY_ORBITAL_PERIOD)
    @test isdefined(ExESS, :MERCURY_PERIHELION)
    @test isdefined(ExESS, :MERCURY_RADIUS)
    @test isdefined(ExESS, :MERCURY_ROTATION_PERIOD)
    @test isdefined(ExESS, :MERCURY_SEMI_MAJOR_AXIS)
    @test isdefined(ExESS, :MERCURY_SYNODIC_ROTATION_PERIOD)

    #::. physical correctness
    @test MERCURY_DAY == MERCURY_SYNODIC_ROTATION_PERIOD
    @test isapprox(MERCURY_SYNODIC_ROTATION_PERIOD, 3 * MERCURY_ROTATION_PERIOD; rtol=1e3)
    @test isapprox(MERCURY_SYNODIC_ROTATION_PERIOD, 2 * MERCURY_ORBITAL_PERIOD; rtol=1e3)
end

println("\rTESTING: mercury > constants.jl - DONE")
nothing
