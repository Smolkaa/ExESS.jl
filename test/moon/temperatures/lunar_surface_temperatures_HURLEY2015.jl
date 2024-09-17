@testset verbose=true "lunar_surface_temperatures_HURLEY2015.jl" begin

    #
    grid = Spherical2DGrid(1, 720, 360)
    @test lunar_surface_temperatures_HURLEY2015(grid) isa Vector{Float64}

end
