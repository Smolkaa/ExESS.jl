@testset verbose=true "grids.jl ................" begin
    print("TESTING: grids > grids.jl")

    ########################################################################################
    #::. FUNCTION TEST - `mapgrid`
    grid_s2d   = Spherical2DGrid(LUNAR_RADIUS, rand(2:45, 2)...)
    grid_s2de  = Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(2:45, 2)...)
    grid_s2dr  = Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(2:45))
    grid_s2dre = Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(2:45))
    x_s2d   = collect(1:length(grid_s2d))
    x_s2de  = collect(1:length(grid_s2de))
    x_s2dr  = collect(1:length(grid_s2dr))
    x_s2dre = collect(1:length(grid_s2dre))

    # test no error
    @test all(_ -> mapgrid(x_s2d, grid_s2d, Spherical2DGrid(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid_s2d, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid_s2d, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid_s2d, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2de, grid_s2de, Spherical2DGrid(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_s2de, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_s2de, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_s2de, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dr, grid_s2dr, Spherical2DGrid(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_s2dr, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_s2dr, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_s2dr, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dre, grid_s2dre, Spherical2DGrid(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_s2dre, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(2:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_s2dre, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_s2dre, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(2:45))) isa Vector{Int64}, 1:100)

    # test invariance if the same grid is used
    @test mapgrid(x_s2d, grid_s2d, grid_s2d) == x_s2d
    @test mapgrid(x_s2de, grid_s2de, grid_s2de) == x_s2de
    @test mapgrid(x_s2dr, grid_s2dr, grid_s2dr) == x_s2dr
    @test mapgrid(x_s2dre, grid_s2dre, grid_s2dre) == x_s2dre

    # test output size
    @test all(g -> length(mapgrid(x_s2d, grid_s2d, g)) == length(g), [grid_s2d, grid_s2de, grid_s2dr, grid_s2dre])
    @test all(g -> length(mapgrid(x_s2de, grid_s2de, g)) == length(g), [grid_s2d, grid_s2de, grid_s2dr, grid_s2dre])
    @test all(g -> length(mapgrid(x_s2dr, grid_s2dr, g)) == length(g), [grid_s2d, grid_s2de, grid_s2dr, grid_s2dre])
    @test all(g -> length(mapgrid(x_s2dre, grid_s2dre, g)) == length(g), [grid_s2d, grid_s2de, grid_s2dr, grid_s2dre])
    ########################################################################################


    println("\rTESTING: grids > grids.jl - DONE")
end

nothing
