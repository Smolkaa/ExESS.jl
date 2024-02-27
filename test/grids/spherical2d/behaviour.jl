@testset verbose=true "behaviour" begin

    # create grids
    r, N_theta, N_phi = rand()*100, rand(2:180), rand(2:90)
    grid = Spherical2DGrid(r, N_theta, N_phi)
    grid_eqsim = Spherical2DGrid_EqSim(r, N_theta, N_phi)
    grid_reduced = Spherical2DGrid_Reduced(r, N_phi)
    grid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(r, N_phi)

    # test lengths
    @test length(grid) == N_theta*N_phi == length(coords(grid)) == length(areas(grid))
    @test length(grid_eqsim) == N_theta*N_phi == length(coords(grid_eqsim)) == length(areas(grid_eqsim))
    @test length(grid_reduced) == length(coords(grid_reduced)) == length(areas(grid_reduced))
    @test length(grid_reduced_eqsim) == length(coords(grid_reduced_eqsim)) == length(areas(grid_reduced_eqsim))

    # test sizes
    @test size(grid) == (N_theta, N_phi)
    @test size(grid_eqsim) == (N_theta, N_phi)
    @test size(grid_reduced)[1] isa AbstractVector
    @test size(grid_reduced_eqsim)[1] isa AbstractVector

    # test completeness
    RTOL = 1e-2
    @test isapprox(sum(areas(grid)), 4*pi*r^2, rtol=RTOL)
    @test isapprox(sum(areas(grid_eqsim)), 2*pi*r^2, rtol=RTOL)
    @test isapprox(sum(areas(grid_reduced)), 4*pi*r^2, rtol=RTOL)
    @test isapprox(sum(areas(grid_reduced_eqsim)), 2*pi*r^2, rtol=RTOL)

    # test for zero volume of 2D grid
    @test sum(volumes(grid)) == 0
    @test sum(volumes(grid_eqsim)) == 0
    @test sum(volumes(grid_reduced)) == 0
    @test sum(volumes(grid_reduced_eqsim)) == 0

    # test for correct surface coordinates (= all coordinates in 2D grids)
    @test surfacecoords(grid) == coords(grid) isa AbstractVector
    @test surfacecoords(grid_eqsim) == coords(grid_eqsim) isa AbstractVector
    @test surfacecoords(grid_reduced) == coords(grid_reduced) isa AbstractVector
    @test surfacecoords(grid_reduced_eqsim) == coords(grid_reduced_eqsim) isa AbstractVector

    # test for correct index and coordinate conversion
    xs_offset = (0, 0.001*rand()*pi, 0.001*rand()*pi/2)
    for _ in 1:100
        idx = rand(1:length(grid))
        xs = coords(grid)[idx]
        @test coord2idx(grid, xs) == idx
        @test coord2idx(grid, xs + xs_offset) == idx

        idx = rand(1:length(grid_eqsim))
        xs = coords(grid_eqsim)[idx]
        @test coord2idx(grid_eqsim, xs) == idx
        @test coord2idx(grid_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(grid_reduced))
        xs = coords(grid_reduced)[idx]
        @test coord2idx(grid_reduced, xs) == idx
        @test coord2idx(grid_reduced, xs + xs_offset) == idx
        
        idx = rand(1:length(grid_reduced_eqsim))
        xs = coords(grid_reduced_eqsim)[idx]
        @test coord2idx(grid_reduced_eqsim, xs) == idx
        @test coord2idx(grid_reduced_eqsim, xs + xs_offset) == idx

    end
end