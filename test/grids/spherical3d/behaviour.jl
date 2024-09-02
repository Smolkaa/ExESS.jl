@testset verbose=true "behaviour" begin

    # create grid_eqsim
    R, N_lon, N_lat = sort(rand(rand(2:10)))*100, rand(2:180), rand(2:90)
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
    @test length(grid) == (length(R)-1)*N_lon*N_lat == length(coords(grid)) == length(areas(grid)) == length(volumes(grid))

    # test sizes
    @test size(grid) == (length(R)-1, N_lon, N_lat)

    # test area completeness
    RTOL = 1e-3
    @test all(i -> isapprox(sum(areas(grid)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), 4*pi*R[i]^2, rtol=RTOL), 1:length(R)-1)

    # test volume completeness
    RTOL = 1e-3
    @test isapprox(sum(volumes(grid)), 4/3*pi*(R[end].^3 - R[1].^3), rtol=RTOL)

    # surface coordinates
    @test surfacecoords(grid) == coords(grid)[1:N_lon*N_lat]

    # test for correct index and coordinate conversion
    xs_offset = (rand(), rand()*pi, rand()*pi/2) .* 1e-5
    for _ in 1:100
        idx = rand(1:length(grid))
        xs = coords(grid)[idx]
        @test coord2idx(grid, xs) == idx
        @test coord2idx(grid, xs + xs_offset) == idx
    end

    # test for overflowing longitude (pclamp to lonrange)
    @test all(x -> coord2idx(grid, x...) == coord2idx(grid, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    # test for overflowing latitude (clamp to max. ranges)
    @test all(x -> coord2idx(grid, x[1], x[2], pi/2 + rand()) == coord2idx(grid, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    # test for NOT drawing indecies < 1 at grid element borders
    @test all(i -> all(j -> all(k -> coord2idx(grid, grid.r0 + grid.h[i], grid.lonrange[j], grid.latrange[k]) > 0, 1:2), 1:2), 1:length(grid.h)-1)

    # test mapgrid - setup
    x_s3d   = collect(1:length(grid))

    # test mapgrid

    # test mapgrid - invariance if the same grid is used
    @test mapgrid(x_s3d, grid, grid) == x_s3d
end
