@testset verbose=true "behaviour" begin

    # create grid_eqsim
    R, N_lon, N_lat = sort(rand(rand(2:10)))*100, rand(2:180), rand(2:90)
    N_R = length(R) - 1
    grid = Spherical3DGrid(R, N_lon, N_lat)
    grid_eqsim = Spherical3DGrid_EqSim(R, N_lon, N_lat)
    grid_reduced = Spherical3DGrid_Reduced(R, N_lat)
    grid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(R, N_lat)

    # create partial grids with random ranges
    lonmin, latmin = rand(2) .* rand([-1,1])
    lonmax, latmax = lonmin + rand()*pi/10, latmin + rand()*pi/10
    pgrid = Spherical3DGrid(R, N_lon, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
    # pgrid_eqsim = Spherical3DGrid_EqSim(R, N_lon, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))
    # pgrid_reduced = Spherical3DGrid_Reduced(R, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
    # pgrid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(R, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))

    # test lengths
    @test length(grid) == N_R*N_lon*N_lat == length(coords(grid)) == length(areas(grid)) == length(volumes(grid))
    @test length(grid_eqsim) == N_R*N_lon*N_lat == length(coords(grid_eqsim)) == length(areas(grid_eqsim)) == length(volumes(grid_eqsim))
    @test length(grid_reduced) == length(coords(grid_reduced)) == length(areas(grid_reduced)) == length(volumes(grid_reduced))

    # test sizes
    @test size(grid) == (N_R, N_lon, N_lat)
    @test size(grid_eqsim) == (N_R, N_lon, N_lat)
    @test size(grid_reduced)[2] isa AbstractVector

    # test area completeness
    RTOL = 1e-3
    Int64(length(grid_reduced) / N_R)
    @test all(i -> isapprox(sum(areas(grid)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), 4*pi*R[i]^2, rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(grid_eqsim)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), 2*pi*R[i]^2, rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(grid_reduced)[1+(i-1)*Int64(length(grid_reduced)/N_R):i*Int64(length(grid_reduced)/N_R)]), 4*pi*R[i]^2, rtol=RTOL), 1:N_R)

    # test volume completeness
    RTOL = 1e-3
    @test isapprox(sum(volumes(grid)), 4/3*pi*(R[end].^3 - R[1].^3), rtol=RTOL)
    @test isapprox(sum(volumes(grid_eqsim)), 2/3*pi*(R[end].^3 - R[1].^3), rtol=RTOL)

    # surface coordinates
    @test surfacecoords(grid) == coords(grid)[1:N_lon*N_lat]
    @test surfacecoords(grid_eqsim) == coords(grid_eqsim)[1:N_lon*N_lat]

    # test for correct index and coordinate conversion
    xs_offset = (rand(), rand()*pi, rand()*pi/2) .* 1e-5
    for _ in 1:100
        idx = rand(1:length(grid))
        xs = coords(grid)[idx]
        @test coord2idx(grid, xs) == idx
        @test coord2idx(grid, xs + xs_offset) == idx

        idx = rand(1:length(grid_eqsim))
        xs = coords(grid_eqsim)[idx]
        @test coord2idx(grid_eqsim, xs) == idx
        @test coord2idx(grid_eqsim, xs + xs_offset) == idx
    end

    # test for overflowing longitude (pclamp to lonrange)
    @test all(x -> coord2idx(grid, x...) == coord2idx(grid, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid, x...) == coord2idx(grid, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    # test for overflowing latitude (clamp to max. ranges)
    @test all(x -> coord2idx(grid, x[1], x[2], pi/2 + rand()) == coord2idx(grid, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid, x[1], x[2], -pi/2 - rand()) == coord2idx(grid, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    # test for NOT drawing indecies < 1 at grid element borders
    @test all(i -> all(j -> all(k -> coord2idx(grid, grid.r0 + grid.h[i], grid.lonrange[j], grid.latrange[k]) > 0, 1:2), 1:2), 1:length(grid.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(grid_eqsim, grid_eqsim.r0 + grid_eqsim.h[i], grid_eqsim.lonrange[j], grid_eqsim.latrange[k]) > 0, 1:2), 1:2), 1:length(grid_eqsim.h)-1)

    # test mapgrid - setup
    x_s3d   = collect(1:length(grid))
    x_s3de  = collect(1:length(grid_eqsim))

    # test mapgrid
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid(grid.r0, grid.h, rand(5:45,2)...)) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_EqSim(grid.r0, grid.h, rand(5:45,2)...)) isa typeof(x_s3d), 1:100)

    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_EqSim(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3de), 1:100)

    # test mapgrid - invariance if the same grid is used
    @test mapgrid(x_s3d, grid, grid) == x_s3d
    @test mapgrid(x_s3de, grid_eqsim, grid_eqsim) == x_s3de
end
