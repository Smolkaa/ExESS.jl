@testset verbose=true "behaviour" begin

    # create grids
    r, N_theta, N_phi = rand()*100, rand(2:180), rand(2:90)
    grid = Spherical2DGrid(r, N_theta, N_phi)
    grid_eqsim = Spherical2DGrid_EqSim(r, N_theta, N_phi)
    grid_reduced = Spherical2DGrid_Reduced(r, N_phi)
    grid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(r, N_phi)

    # test lengths and sizes
    @test length(grid) == N_theta*N_phi == length(coords(grid)) == length(areas(grid))
    @test length(grid_eqsim) == N_theta*N_phi == length(coords(grid_eqsim)) == length(areas(grid_eqsim))
    @test length(grid_reduced) == length(coords(grid_reduced)) == length(areas(grid_reduced))
    @test length(grid_reduced_eqsim) == length(coords(grid_reduced_eqsim)) == length(areas(grid_reduced_eqsim))

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

    # TODO: surfacecoords, coord2idx,...
end