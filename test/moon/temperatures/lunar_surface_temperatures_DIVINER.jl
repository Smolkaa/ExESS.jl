@testset verbose=true "lunar_surface_temperatures_DIVINER.jl" begin

    #
    grid = Spherical2DGrid(1, 720, 360)
    @test lunar_surface_temperatures_DIVINER(rand(), grid) isa Vector{Float64}
    @test lunar_surface_temperatures_DIVINER_avg( grid) isa Vector{Float64}

end
