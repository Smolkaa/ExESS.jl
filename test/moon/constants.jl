@testset verbose=true "constants.jl .................." begin
    @test isdefined(ExESS, :LUNAR_RADIUS)
    @test isdefined(ExESS, :LUNAR_MASS)
    @test isdefined(ExESS, :LUNAR_DAY)
    @test isdefined(ExESS, :LUNAR_g0)
end