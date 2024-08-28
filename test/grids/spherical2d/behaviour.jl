@testset verbose=true "behaviour" begin

    # create grids
    r, N_theta, N_phi = rand()*100, rand(2:180), rand(2:90)
    grid = Spherical2DGrid(r, N_theta, N_phi)
    grid_eqsim = Spherical2DGrid_EqSim(r, N_theta, N_phi)
    grid_reduced = Spherical2DGrid_Reduced(r, N_phi)
    grid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(r, N_phi)

    # create partial grids with random ranges
    lonmin, latmin = rand(2) .* rand([-1,1])
    lonmax, latmax = lonmin + rand()*pi/10, latmin + rand()*pi/10
    pgrid = Spherical2DGrid(r, N_theta, N_phi; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
    pgrid_eqsim = Spherical2DGrid_EqSim(r, N_theta, N_phi; lonrange=(lonmin, lonmax), latmax=abs(latmax))
    pgrid_reduced = Spherical2DGrid_Reduced(r, N_phi; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
    pgrid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(r, N_phi; lonrange=(lonmin, lonmax), latmax=abs(latmax))

    # test lengths
    @test length(grid) == N_theta*N_phi == length(coords(grid)) == length(areas(grid))
    @test length(grid_eqsim) == N_theta*N_phi == length(coords(grid_eqsim)) == length(areas(grid_eqsim))
    @test length(grid_reduced) == length(coords(grid_reduced)) == length(areas(grid_reduced))
    @test length(grid_reduced_eqsim) == length(coords(grid_reduced_eqsim)) == length(areas(grid_reduced_eqsim))

    @test length(pgrid) == N_theta*N_phi == length(coords(pgrid)) == length(areas(pgrid))
    @test length(pgrid_eqsim) == N_theta*N_phi == length(coords(pgrid_eqsim)) == length(areas(pgrid_eqsim))
    @test length(pgrid_reduced) == length(coords(pgrid_reduced)) == length(areas(pgrid_reduced))
    @test length(pgrid_reduced_eqsim) == length(coords(pgrid_reduced_eqsim)) == length(areas(pgrid_reduced_eqsim))

    # test sizes
    @test size(grid) == (N_theta, N_phi)
    @test size(grid_eqsim) == (N_theta, N_phi)
    @test size(grid_reduced)[1] isa AbstractVector
    @test size(grid_reduced_eqsim)[1] isa AbstractVector

    @test size(pgrid) == (N_theta, N_phi)
    @test size(pgrid_eqsim) == (N_theta, N_phi)
    @test size(pgrid_reduced)[1] isa AbstractVector
    @test size(pgrid_reduced_eqsim)[1] isa AbstractVector

    # test completeness
    RTOL = 1e-2
    @test isapprox(sum(areas(grid)), 4*pi*r^2; rtol=RTOL)
    @test isapprox(sum(areas(grid_eqsim)), 2*pi*r^2; rtol=RTOL)
    @test isapprox(sum(areas(grid_reduced)), 4*pi*r^2; rtol=RTOL)
    @test isapprox(sum(areas(grid_reduced_eqsim)), 2*pi*r^2; rtol=RTOL)

    @test isapprox(sum(areas(pgrid)), r^2 * (lonmax - lonmin) * (sin(latmax) - sin(latmin)); rtol=RTOL)
    @test isapprox(sum(areas(pgrid_eqsim)), r^2 * (lonmax - lonmin) * sin(abs(latmax)); rtol=RTOL)
    @test isapprox(sum(areas(pgrid_reduced)), r^2 * (lonmax - lonmin) * (sin(latmax) - sin(latmin)); rtol=RTOL)
    @test isapprox(sum(areas(pgrid_reduced_eqsim)), r^2 * (lonmax - lonmin) * sin(abs(latmax)); rtol=RTOL)

    # test for zero volume of 2D grid
    @test sum(volumes(grid)) == 0
    @test sum(volumes(grid_eqsim)) == 0
    @test sum(volumes(grid_reduced)) == 0
    @test sum(volumes(grid_reduced_eqsim)) == 0

    @test sum(volumes(pgrid)) == 0
    @test sum(volumes(pgrid_eqsim)) == 0
    @test sum(volumes(pgrid_reduced)) == 0
    @test sum(volumes(pgrid_reduced_eqsim)) == 0

    # test for correct surface coordinates (= all coordinates in 2D grids)
    @test surfacecoords(grid) == coords(grid) isa AbstractVector
    @test surfacecoords(grid_eqsim) == coords(grid_eqsim) isa AbstractVector
    @test surfacecoords(grid_reduced) == coords(grid_reduced) isa AbstractVector
    @test surfacecoords(grid_reduced_eqsim) == coords(grid_reduced_eqsim) isa AbstractVector

    @test surfacecoords(pgrid) == coords(pgrid) isa AbstractVector
    @test surfacecoords(pgrid_eqsim) == coords(pgrid_eqsim) isa AbstractVector
    @test surfacecoords(pgrid_reduced) == coords(pgrid_reduced) isa AbstractVector
    @test surfacecoords(pgrid_reduced_eqsim) == coords(pgrid_reduced_eqsim) isa AbstractVector

    # test for correct index and coordinate conversion
    xs_offset = (0, 0.0001*rand()*pi, 0.0001*rand()*pi/2)
    for _ in 1:100
        idx = rand(1:length(grid))
        xs = coords(grid)[idx]
        @test coord2idx(grid, xs) == idx
        @test coord2idx(grid, xs + xs_offset) == idx

        idx = rand(1:length(grid_eqsim))
        xs = coords(grid_eqsim)[idx]
        @test coord2idx(grid_eqsim, xs) == idx
        @test coord2idx(grid_eqsim, xs.theta, -xs.phi) == idx
        @test coord2idx(grid_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(grid_reduced))
        xs = coords(grid_reduced)[idx]
        @test coord2idx(grid_reduced, xs) == idx
        @test coord2idx(grid_reduced, xs + xs_offset) == idx

        idx = rand(1:length(grid_reduced_eqsim))
        xs = coords(grid_reduced_eqsim)[idx]
        @test coord2idx(grid_reduced_eqsim, xs) == idx
        @test coord2idx(grid_reduced_eqsim, xs.theta, -xs.phi) == idx
        @test coord2idx(grid_reduced_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(pgrid))
        xs = coords(pgrid)[idx]
        @test coord2idx(pgrid, xs) == idx
        @test coord2idx(pgrid, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_eqsim))
        xs = coords(pgrid_eqsim)[idx]
        @test coord2idx(pgrid_eqsim, xs) == idx
        @test coord2idx(pgrid_eqsim, xs.theta, -xs.phi) == idx
        @test coord2idx(pgrid_eqsim, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_reduced))
        xs = coords(pgrid_reduced)[idx]
        @test coord2idx(pgrid_reduced, xs) == idx
        @test coord2idx(pgrid_reduced, xs + xs_offset) == idx

        idx = rand(1:length(pgrid_reduced_eqsim))
        xs = coords(pgrid_reduced_eqsim)[idx]
        @test coord2idx(pgrid_reduced_eqsim, xs) == idx
        @test coord2idx(pgrid_reduced_eqsim, xs.theta, -xs.phi) == idx
        @test coord2idx(pgrid_reduced_eqsim, xs + xs_offset) == idx
    end

    # test for overflowing longitudes (periodic clamping to [-pi, pi])
    @test all(r -> coord2idx(grid, r...) == coord2idx(grid, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_eqsim, r...) == coord2idx(grid_eqsim, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_reduced, r...) == coord2idx(grid_reduced, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_reduced_eqsim, r...) == coord2idx(grid_reduced_eqsim, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])

    @test all(r -> coord2idx(pgrid, r...) == coord2idx(pgrid, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_eqsim, r...) == coord2idx(pgrid_eqsim, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_reduced, r...) == coord2idx(pgrid_reduced, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_reduced_eqsim, r...) == coord2idx(pgrid_reduced_eqsim, 2pi + r[1], r[2]), [rand(2) for _ in 1:1000])

    @test all(r -> coord2idx(grid, -r[1], r[2]) == coord2idx(grid, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_eqsim, -r[1], r[2]) == coord2idx(grid_eqsim, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_reduced, -r[1], r[2]) == coord2idx(grid_reduced, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(grid_reduced_eqsim, -r[1], r[2]) == coord2idx(grid_reduced_eqsim, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])

    @test all(r -> coord2idx(pgrid, -r[1], r[2]) == coord2idx(pgrid, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_eqsim, -r[1], r[2]) == coord2idx(pgrid_eqsim, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_reduced, -r[1], r[2]) == coord2idx(pgrid_reduced, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])
    @test all(r -> coord2idx(pgrid_reduced_eqsim, -r[1], r[2]) == coord2idx(pgrid_reduced_eqsim, -2pi - r[1], r[2]), [rand(2) for _ in 1:1000])

    # test coord2idx for overflowing latitudes (clamping to max. ranges)
    @test all(r -> coord2idx(grid, r, pi/2 + rand()) == coord2idx(grid, r, pi/2), rand(1000))
    @test all(r -> coord2idx(grid_eqsim, r, pi/2 + rand()) == coord2idx(grid_eqsim, r, pi/2), rand(1000))
    @test all(r -> coord2idx(grid_reduced, r, pi/2 + rand()) == coord2idx(grid_reduced, r, pi/2), rand(1000))
    @test all(r -> coord2idx(grid_reduced_eqsim, r, pi/2 + rand()) == coord2idx(grid_reduced_eqsim, r, pi/2), rand(1000))

    @test all(r -> coord2idx(pgrid, r, pi/2 + rand()) == coord2idx(pgrid, r, pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_eqsim, r, pi/2 + rand()) == coord2idx(pgrid_eqsim, r, pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_reduced, r, pi/2 + rand()) == coord2idx(pgrid_reduced, r, pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_reduced_eqsim, r, pi/2 + rand()) == coord2idx(pgrid_reduced_eqsim, r, pi/2), rand(1000))

    @test all(r -> coord2idx(grid, r, -pi/2 - rand()) == coord2idx(grid, r, -pi/2), rand(1000))
    @test all(r -> coord2idx(grid_eqsim, r, -pi/2 - rand()) == coord2idx(grid_eqsim, r, pi/2), rand(1000))
    @test all(r -> coord2idx(grid_reduced, r, -pi/2 - rand()) == coord2idx(grid_reduced, r, -pi/2), rand(1000))
    @test all(r -> coord2idx(grid_reduced_eqsim, r, -pi/2 - rand()) == coord2idx(grid_reduced_eqsim, r, pi/2), rand(1000))

    @test all(r -> coord2idx(pgrid, r, -pi/2 - rand()) == coord2idx(pgrid, r, -pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_eqsim, r, -pi/2 - rand()) == coord2idx(pgrid_eqsim, r, pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_reduced, r, -pi/2 - rand()) == coord2idx(pgrid_reduced, r, -pi/2), rand(1000))
    @test all(r -> coord2idx(pgrid_reduced_eqsim, r, -pi/2 - rand()) == coord2idx(pgrid_reduced_eqsim, r, pi/2), rand(1000))

    # test for NOT drawing index <1 at grid element borders
    # NOTE: using `BigFloat` does not work here, though it should never be used
    EPS = [eps(Float16), eps(Float32), eps(Float64), -eps(Float16), -eps(Float32), -eps(Float64)]
    for g in [grid, grid_eqsim, grid_reduced, grid_reduced_eqsim]
        for tshift in EPS
            for pshift in EPS
                @test coord2idx(g, g.lonrange[1] + tshift, g.latrange[1] + pshift) > 0
                @test coord2idx(g, g.lonrange[2] + tshift, g.latrange[1] + pshift) > 0
                @test coord2idx(g, g.lonrange[1] + tshift, g.latrange[2] + pshift) > 0
                @test coord2idx(g, g.lonrange[2] + tshift, g.latrange[2] + pshift) > 0
            end
        end
        # @test all(tshift -> all(pshift -> coord2idx(g, g.lonrange[1] + tshift, g.latrange[1] + pshift) > 0, EPS), EPS)
        # @test all(tshift -> all(pshift -> coord2idx(g, g.lonrange[2] + tshift, g.latrange[1] + pshift) > 0, EPS), EPS)
        # @test all(tshift -> all(pshift -> coord2idx(g, g.lonrange[1] + tshift, g.latrange[2] + pshift) > 0, EPS), EPS)
        # @test all(tshift -> all(pshift -> coord2idx(g, g.lonrange[2] + tshift, g.latrange[2] + pshift) > 0, EPS), EPS)
    end

    # test mapgrid - setup
    x_s2d   = collect(1:length(grid))
    x_s2de  = collect(1:length(grid_eqsim))
    x_s2dr  = collect(1:length(grid_reduced))
    x_s2dre = collect(1:length(grid_reduced_eqsim))

    x_ps2d   = collect(1:length(pgrid))
    x_ps2de  = collect(1:length(pgrid_eqsim))
    x_ps2dr  = collect(1:length(pgrid_reduced))
    x_ps2dre = collect(1:length(pgrid_reduced_eqsim))

    # test mapgrid - behaviour
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2d, grid, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2de, grid_eqsim, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dr, grid_reduced, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45))) isa Vector{Int64}, 1:100)

    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_EqSim(LUNAR_RADIUS, rand(5:45, 2)...;
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_Reduced(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latrange=(0.5*rand(), 0.5 + rand()*0.5))) isa Vector{Int64}, 1:100)
    @test all(_ -> mapgrid(x_s2dre, grid_reduced_eqsim, Spherical2DGrid_Reduced_EqSim(LUNAR_RADIUS, rand(5:45);
        lonrange=(rand(), 1+rand()), latmax=0.5 + rand()*0.5)) isa Vector{Int64}, 1:100)

    @test_throws BoundsError mapgrid(x_ps2d, pgrid, grid)
    @test_throws BoundsError mapgrid(x_ps2d, pgrid, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps2d, pgrid, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps2d, pgrid, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps2de, pgrid_eqsim, grid)
    @test_throws BoundsError mapgrid(x_ps2de, pgrid_eqsim, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps2de, pgrid_eqsim, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps2de, pgrid_eqsim, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps2dr, pgrid_reduced, grid)
    @test_throws BoundsError mapgrid(x_ps2dr, pgrid_reduced, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps2dr, pgrid_reduced, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps2dr, pgrid_reduced, grid_reduced_eqsim)

    @test_throws BoundsError mapgrid(x_ps2dre, pgrid_reduced_eqsim, grid)
    @test_throws BoundsError mapgrid(x_ps2dre, pgrid_reduced_eqsim, grid_eqsim)
    @test_throws BoundsError mapgrid(x_ps2dre, pgrid_reduced_eqsim, grid_reduced)
    @test_throws BoundsError mapgrid(x_ps2dre, pgrid_reduced_eqsim, grid_reduced_eqsim)


    # test mapgrid - invariance if the same grid is used
    @test mapgrid(x_s2d, grid, grid) == x_s2d
    @test mapgrid(x_s2de, grid_eqsim, grid_eqsim) == x_s2de
    @test mapgrid(x_s2dr, grid_reduced, grid_reduced) == x_s2dr
    @test mapgrid(x_s2dre, grid_reduced_eqsim, grid_reduced_eqsim) == x_s2dre

    @test mapgrid(x_ps2d, pgrid, pgrid) == x_ps2d
    @test mapgrid(x_ps2de, pgrid_eqsim, pgrid_eqsim) == x_ps2de
    @test mapgrid(x_ps2dr, pgrid_reduced, pgrid_reduced) == x_ps2dr
    @test mapgrid(x_ps2dre, pgrid_reduced_eqsim, pgrid_reduced_eqsim) == x_ps2dre

    # test mapgrid - output size
    G = [grid, grid_eqsim, grid_reduced, grid_reduced_eqsim,
         pgrid, pgrid_eqsim, pgrid_reduced, pgrid_reduced_eqsim]
    @test all(g -> length(mapgrid(x_s2d, grid, g)) == length(g), G)
    @test all(g -> length(mapgrid(x_s2de, grid_eqsim, g)) == length(g), G)
    @test all(g -> length(mapgrid(x_s2dr, grid_reduced, g)) == length(g), G)
    @test all(g -> length(mapgrid(x_s2dre, grid_reduced_eqsim, g)) == length(g), G)
end
