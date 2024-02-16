@testset verbose=true "constants.jl .................." begin
    @test isdefined(ExESS, :AVOGADRO_CONSTANT)
    @test isdefined(ExESS, :BOLTZMANN_CONSTANT)
    @test isdefined(ExESS, :GRAVITATIONAL_CONSTANT)
    @test isdefined(ExESS, :PLANCK_CONSTANT)
    @test isdefined(ExESS, :UNIVERSAL_GAS_CONSTANT)
    @test isdefined(ExESS, :ELEMENTARY_CHARGE)

    @test isdefined(ExESS, :AMU_H)
    @test isdefined(ExESS, :AMU_He)
    @test isdefined(ExESS, :AMU_O)
    @test isdefined(ExESS, :AMU_Ne)
    @test isdefined(ExESS, :AMU_Ar)

    @test isdefined(ExESS, :AMU_H2)
    @test isdefined(ExESS, :AMU_OH)
    @test isdefined(ExESS, :AMU_H2O)

    @test isdefined(ExESS, :MERCURY_MASS)
    @test isdefined(ExESS, :MERCURY_RADIUS)

    @test isdefined(ExESS, :CERES_MASS)
    @test isdefined(ExESS, :CERES_RADIUS)
end