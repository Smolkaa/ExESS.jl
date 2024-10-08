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
    pgrid_eqsim = Spherical3DGrid_EqSim(R, N_lon, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))
    pgrid_reduced = Spherical3DGrid_Reduced(R, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
    pgrid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(R, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))


    # test lengths
    @test length(grid) == N_R*N_lon*N_lat == length(coords(grid)) == length(areas(grid)) == length(volumes(grid))
    @test length(grid_eqsim) == N_R*N_lon*N_lat == length(coords(grid_eqsim)) == length(areas(grid_eqsim)) == length(volumes(grid_eqsim))
    @test length(grid_reduced) == length(coords(grid_reduced)) == length(areas(grid_reduced)) == length(volumes(grid_reduced))
    @test length(grid_reduced_eqsim) == length(coords(grid_reduced_eqsim)) == length(areas(grid_reduced_eqsim)) == length(volumes(grid_reduced_eqsim))

    @test length(pgrid) == N_R*N_lon*N_lat == length(coords(pgrid)) == length(areas(pgrid)) == length(volumes(pgrid))
    @test length(pgrid_eqsim) == N_R*N_lon*N_lat == length(coords(pgrid_eqsim)) == length(areas(pgrid_eqsim)) == length(volumes(pgrid_eqsim))
    @test length(pgrid_reduced) == length(coords(pgrid_reduced)) == length(areas(pgrid_reduced)) == length(volumes(pgrid_reduced))
    @test length(pgrid_reduced_eqsim) == length(coords(pgrid_reduced_eqsim)) == length(areas(pgrid_reduced_eqsim)) == length(volumes(pgrid_reduced_eqsim))


    # test sizes
    @test size(grid) == (N_R, N_lon, N_lat)
    @test size(grid_eqsim) == (N_R, N_lon, N_lat)
    @test size(grid_reduced)[2] isa AbstractVector
    @test size(grid_reduced_eqsim)[2] isa AbstractVector

    @test size(pgrid) == (N_R, N_lon, N_lat)
    @test size(pgrid_eqsim) == (N_R, N_lon, N_lat)
    @test size(pgrid_reduced)[2] isa AbstractVector
    @test size(pgrid_reduced_eqsim)[2] isa AbstractVector


    # test area completeness
    RTOL = 1e-4
    Int64(length(grid_reduced) / N_R)
    @test all(i -> isapprox(sum(areas(grid)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), 4*pi*R[i]^2, rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(grid_eqsim)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), 2*pi*R[i]^2, rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(grid_reduced)[1+(i-1)*Int64(length(grid_reduced)/N_R):i*Int64(length(grid_reduced)/N_R)]), 4*pi*R[i]^2, rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(grid_reduced_eqsim)[1+(i-1)*Int64(length(grid_reduced_eqsim)/N_R):i*Int64(length(grid_reduced_eqsim)/N_R)]), 2*pi*R[i]^2, rtol=RTOL), 1:N_R)

    @test all(i -> isapprox(sum(areas(pgrid)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), R[i]^2*(lonmax-lonmin)*(sin(latmax)-sin(latmin)), rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(pgrid_eqsim)[1+(i-1)*N_lon*N_lat:i*N_lon*N_lat]), R[i]^2*(lonmax-lonmin)*sin(abs(latmax)), rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(pgrid_reduced)[1+(i-1)*Int64(length(pgrid_reduced)/N_R):i*Int64(length(pgrid_reduced)/N_R)]), R[i]^2*(lonmax-lonmin)*(sin(latmax)-sin(latmin)), rtol=RTOL), 1:N_R)
    @test all(i -> isapprox(sum(areas(pgrid_reduced_eqsim)[1+(i-1)*Int64(length(pgrid_reduced_eqsim)/N_R):i*Int64(length(pgrid_reduced_eqsim)/N_R)]), R[i]^2*(lonmax-lonmin)*sin(abs(latmax)), rtol=RTOL), 1:N_R)


    # test volume completeness
    RTOL = 1e-3
    @test isapprox(sum(volumes(grid)), 4/3*pi*(R[end]^3 - R[1]^3), rtol=RTOL)
    @test isapprox(sum(volumes(grid_eqsim)), 2/3*pi*(R[end]^3 - R[1]^3), rtol=RTOL)
    @test isapprox(sum(volumes(grid_reduced)), 4/3*pi*(R[end]^3 - R[1]^3), rtol=RTOL)
    @test isapprox(sum(volumes(grid_reduced_eqsim)), 2/3*pi*(R[end]^3 - R[1]^3), rtol=RTOL)

    @test isapprox(sum(volumes(pgrid)), (R[end]^3 - R[1]^3)/3 * (lonmax - lonmin) * (sin(latmax) - sin(latmin)), rtol=RTOL)
    @test isapprox(sum(volumes(pgrid_eqsim)), (R[end]^3 - R[1]^3)/3 * (lonmax - lonmin) * sin(abs(latmax)), rtol=RTOL)
    @test isapprox(sum(volumes(pgrid_reduced)), (R[end]^3 - R[1]^3)/3 * (lonmax - lonmin) * (sin(latmax) - sin(latmin)), rtol=RTOL)
    @test isapprox(sum(volumes(pgrid_reduced_eqsim)), (R[end]^3 - R[1]^3)/3 * (lonmax - lonmin) * sin(abs(latmax)), rtol=RTOL)


    # surface coordinates
    @test surfacecoords(grid) == coords(grid)[1:N_lon*N_lat]
    @test surfacecoords(grid_eqsim) == coords(grid_eqsim)[1:N_lon*N_lat]
    @test surfacecoords(grid_reduced) == coords(grid_reduced)[1:Int64(length(grid_reduced)/N_R)]
    @test surfacecoords(grid_reduced_eqsim) == coords(grid_reduced_eqsim)[1:Int64(length(grid_reduced_eqsim)/N_R)]

    @test surfacecoords(pgrid) == coords(pgrid)[1:N_lon*N_lat]
    @test surfacecoords(pgrid_eqsim) == coords(pgrid_eqsim)[1:N_lon*N_lat]
    @test surfacecoords(pgrid_reduced) == coords(pgrid_reduced)[1:Int64(length(pgrid_reduced)/N_R)]
    @test surfacecoords(pgrid_reduced_eqsim) == coords(pgrid_reduced_eqsim)[1:Int64(length(pgrid_reduced_eqsim)/N_R)]


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
        @test coord2idx(grid_eqsim, xs.r, xs.theta, -xs.phi) == idx
        @test coord2idx(grid_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(grid_reduced))
        xs = coords(grid_reduced)[idx]
        @test coord2idx(grid_reduced, xs) == idx
        @test coord2idx(grid_reduced, xs + xs_offset) == idx

        idx = rand(1:length(grid_reduced_eqsim))
        xs = coords(grid_reduced_eqsim)[idx]
        @test coord2idx(grid_reduced_eqsim, xs) == idx
        @test coord2idx(grid_reduced_eqsim, xs.r, xs.theta, -xs.phi) == idx
        @test coord2idx(grid_reduced_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(pgrid))
        xs = coords(pgrid)[idx]
        @test coord2idx(pgrid, xs) == idx
        @test coord2idx(pgrid, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_eqsim))
        xs = coords(pgrid_eqsim)[idx]
        @test coord2idx(pgrid_eqsim, xs) == idx
        @test coord2idx(pgrid_eqsim, xs.r, xs.theta, -xs.phi) == idx
        @test coord2idx(pgrid_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_reduced))
        xs = coords(pgrid_reduced)[idx]
        @test coord2idx(pgrid_reduced, xs) == idx
        @test coord2idx(pgrid_reduced, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_reduced_eqsim))
        xs = coords(pgrid_reduced_eqsim)[idx]
        @test coord2idx(pgrid_reduced_eqsim, xs) == idx
        @test coord2idx(pgrid_reduced_eqsim, xs.r, xs.theta, -xs.phi) == idx
        @test coord2idx(pgrid_reduced_eqsim, xs + xs_offset) == idx
    end


    # test for overflowing longitude (pclamp to lonrange)
    @test all(x -> coord2idx(grid, x...) == coord2idx(grid, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_eqsim, x...) == coord2idx(grid_eqsim, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced, x...) == coord2idx(grid_reduced, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced_eqsim, x...) == coord2idx(grid_reduced_eqsim, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(pgrid, x...) == coord2idx(pgrid, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_eqsim, x...) == coord2idx(pgrid_eqsim, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced, x...) == coord2idx(pgrid_reduced, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced_eqsim, x...) == coord2idx(pgrid_reduced_eqsim, x[1], x[2] + 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(grid, x...) == coord2idx(grid, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_eqsim, x...) == coord2idx(grid_eqsim, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced, x...) == coord2idx(grid_reduced, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced_eqsim, x...) == coord2idx(grid_reduced_eqsim, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(pgrid, x...) == coord2idx(pgrid, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_eqsim, x...) == coord2idx(pgrid_eqsim, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced, x...) == coord2idx(pgrid_reduced, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced_eqsim, x...) == coord2idx(pgrid_reduced_eqsim, x[1], x[2] - 2pi, x[3]), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])


    # test for overflowing latitude (clamp to max. ranges)
    @test all(x -> coord2idx(grid, x[1], x[2], pi/2 + rand()) == coord2idx(grid, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_eqsim, x[1], x[2], pi/2 + rand()) == coord2idx(grid_eqsim, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced, x[1], x[2], pi/2 + rand()) == coord2idx(grid_reduced, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced_eqsim, x[1], x[2], pi/2 + rand()) == coord2idx(grid_reduced_eqsim, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(pgrid, x[1], x[2], pi/2 + rand()) == coord2idx(pgrid, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_eqsim, x[1], x[2], pi/2 + rand()) == coord2idx(pgrid_eqsim, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced, x[1], x[2], pi/2 + rand()) == coord2idx(pgrid_reduced, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced_eqsim, x[1], x[2], pi/2 + rand()) == coord2idx(pgrid_reduced_eqsim, x[1], x[2], pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(grid, x[1], x[2], -pi/2 - rand()) == coord2idx(grid, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_eqsim, x[1], x[2], -pi/2 - rand()) == coord2idx(grid_eqsim, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced, x[1], x[2], -pi/2 - rand()) == coord2idx(grid_reduced, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(grid_reduced_eqsim, x[1], x[2], -pi/2 - rand()) == coord2idx(grid_reduced_eqsim, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])

    @test all(x -> coord2idx(pgrid, x[1], x[2], -pi/2 - rand()) == coord2idx(pgrid, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_eqsim, x[1], x[2], -pi/2 - rand()) == coord2idx(pgrid_eqsim, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced, x[1], x[2], -pi/2 - rand()) == coord2idx(pgrid_reduced, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])
    @test all(x -> coord2idx(pgrid_reduced_eqsim, x[1], x[2], -pi/2 - rand()) == coord2idx(pgrid_reduced_eqsim, x[1], x[2], -pi/2), [rand(3) .* [R[end]-R[1], 1, 1] .+ [R[1], 0, 0] for _ in 1:1000])


    # test for NOT drawing indecies < 1 at grid grid borders
    @test all(i -> all(j -> all(k -> coord2idx(grid, grid.r0 + grid.h[i], grid.lonrange[j], grid.latrange[k]) > 0, 1:2), 1:2), 1:length(grid.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(grid_eqsim, grid_eqsim.r0 + grid_eqsim.h[i], grid_eqsim.lonrange[j], grid_eqsim.latrange[k]) > 0, 1:2), 1:2), 1:length(grid_eqsim.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(grid_reduced, grid_reduced.r0 + grid_reduced.h[i], grid_reduced.lonrange[j], grid_reduced.latrange[k]) > 0, 1:2), 1:2), 1:length(grid_reduced.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(grid_reduced_eqsim, grid_reduced_eqsim.r0 + grid_reduced_eqsim.h[i], grid_reduced_eqsim.lonrange[j], grid_reduced_eqsim.latrange[k]) > 0, 1:2), 1:2), 1:length(grid_reduced_eqsim.h)-1)

    @test all(i -> all(j -> all(k -> coord2idx(pgrid, pgrid.r0 + pgrid.h[i], pgrid.lonrange[j], pgrid.latrange[k]) > 0, 1:2), 1:2), 1:length(pgrid.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(pgrid_eqsim, pgrid_eqsim.r0 + pgrid_eqsim.h[i], pgrid_eqsim.lonrange[j], pgrid_eqsim.latrange[k]) > 0, 1:2), 1:2), 1:length(pgrid_eqsim.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(pgrid_reduced, pgrid_reduced.r0 + pgrid_reduced.h[i], pgrid_reduced.lonrange[j], pgrid_reduced.latrange[k]) > 0, 1:2), 1:2), 1:length(pgrid_reduced.h)-1)
    @test all(i -> all(j -> all(k -> coord2idx(pgrid_reduced_eqsim, pgrid_reduced_eqsim.r0 + pgrid_reduced_eqsim.h[i], pgrid_reduced_eqsim.lonrange[j], pgrid_reduced_eqsim.latrange[k]) > 0, 1:2), 1:2), 1:length(pgrid_reduced_eqsim.h)-1)


    # test mapgrid - setup
    x_s3d   = collect(1:length(grid))
    x_s3de  = collect(1:length(grid_eqsim))
    x_s3dr  = collect(1:length(grid_reduced))
    x_s3dre = collect(1:length(grid_reduced_eqsim))

    x_ps3d   = collect(1:length(pgrid))
    x_ps3de  = collect(1:length(pgrid_eqsim))
    x_ps3dr  = collect(1:length(pgrid_reduced))
    x_ps3dre = collect(1:length(pgrid_reduced_eqsim))


    # test mapgrid
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid(grid.r0, grid.h, rand(5:45,2)...)) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_EqSim(grid.r0, grid.h, rand(5:45,2)...)) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_Reduced(grid.r0, grid.h, rand(5:45))) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_Reduced_EqSim(grid.r0, grid.h, rand(5:45))) isa typeof(x_s3d), 1:100)

    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid(grid.r0, grid.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_EqSim(grid.r0, grid.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_Reduced(grid.r0, grid.h, rand(5:45); lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3d), 1:100)
    @test all(_ -> mapgrid(x_s3d, grid, Spherical3DGrid_Reduced_EqSim(grid.r0, grid.h, rand(5:45); lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3d), 1:100)

    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_EqSim(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_Reduced(grid_eqsim.r0, grid_eqsim.h, rand(5:45))) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_Reduced_EqSim(grid_eqsim.r0, grid_eqsim.h, rand(5:45))) isa typeof(x_s3de), 1:100)

    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_EqSim(grid_eqsim.r0, grid_eqsim.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_Reduced(grid_eqsim.r0, grid_eqsim.h, rand(5:45); lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3de), 1:100)
    @test all(_ -> mapgrid(x_s3de, grid_eqsim, Spherical3DGrid_Reduced_EqSim(grid_eqsim.r0, grid_eqsim.h, rand(5:45); lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3de), 1:100)

    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid(grid_reduced.r0, grid_reduced.h, rand(5:45,2)...)) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_EqSim(grid_reduced.r0, grid_reduced.h, rand(5:45,2)...)) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_Reduced(grid_reduced.r0, grid_reduced.h, rand(5:45))) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_Reduced_EqSim(grid_reduced.r0, grid_reduced.h, rand(5:45))) isa typeof(x_s3dr), 1:100)

    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid(grid_reduced.r0, grid_reduced.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_EqSim(grid_reduced.r0, grid_reduced.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_Reduced(grid_reduced.r0, grid_reduced.h, rand(5:45); lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3dr), 1:100)
    @test all(_ -> mapgrid(x_s3dr, grid_reduced, Spherical3DGrid_Reduced_EqSim(grid_reduced.r0, grid_reduced.h, rand(5:45); lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3dr), 1:100)

    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_EqSim(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45,2)...)) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_Reduced(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45))) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_Reduced_EqSim(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45))) isa typeof(x_s3dre), 1:100)

    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_EqSim(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45,2)...; lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_Reduced(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45); lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa typeof(x_s3dre), 1:100)
    @test all(_ -> mapgrid(x_s3dre, grid_reduced_eqsim, Spherical3DGrid_Reduced_EqSim(grid_reduced_eqsim.r0, grid_reduced_eqsim.h, rand(5:45); lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa typeof(x_s3dre), 1:100)

    @test_throws BoundsError mapgrid(x_ps3d, pgrid, grid)
    @test_throws BoundsError mapgrid(x_ps3d, pgrid, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps3d, pgrid, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps3d, pgrid, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps3de, pgrid_eqsim, grid)
    @test_throws BoundsError mapgrid(x_ps3de, pgrid_eqsim, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps3de, pgrid_eqsim, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps3de, pgrid_eqsim, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps3dr, pgrid_reduced, grid)
    @test_throws BoundsError mapgrid(x_ps3dr, pgrid_reduced, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps3dr, pgrid_reduced, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps3dr, pgrid_reduced, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps3dre, pgrid_reduced_eqsim, grid)
    @test_throws BoundsError mapgrid(x_ps3dre, pgrid_reduced_eqsim, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps3dre, pgrid_reduced_eqsim, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps3dre, pgrid_reduced_eqsim, grid_reduced_eqsim)


    # test mapgrid - invariance if the same grid is used
    @test mapgrid(x_s3d, grid, grid) == x_s3d
    @test mapgrid(x_s3de, grid_eqsim, grid_eqsim) == x_s3de
    @test mapgrid(x_s3dr, grid_reduced, grid_reduced) == x_s3dr
    @test mapgrid(x_s3dre, grid_reduced_eqsim, grid_reduced_eqsim) == x_s3dre

    @test mapgrid(x_ps3d, pgrid, pgrid) == x_ps3d
    @test mapgrid(x_ps3de, pgrid_eqsim, pgrid_eqsim) == x_ps3de
    @test mapgrid(x_ps3dr, pgrid_reduced, pgrid_reduced) == x_ps3dr
    @test mapgrid(x_ps3dre, pgrid_reduced_eqsim, pgrid_reduced_eqsim) == x_ps3dre


end
