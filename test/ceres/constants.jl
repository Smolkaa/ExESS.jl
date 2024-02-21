@testset verbose=true "constants.jl .................." begin
    @test isdefined(ExESS, :CERES_MASS)
    @test isdefined(ExESS, :CERES_RADIUS)
end