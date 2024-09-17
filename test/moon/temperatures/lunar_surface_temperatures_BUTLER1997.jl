@testset verbose=true "lunar_surface_temperatures_BUTLER1997.jl" begin

    #
    grid = Spherical2DGrid(1, 720, 360)
    @test lunar_surface_temperatures_BUTLER1997(grid) isa Vector{Float64}

end
