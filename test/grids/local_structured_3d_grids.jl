# TODO: REWORK TESTS >> GRID TYPES OUTSIDE OF TYPE-LOOP
# TODO: HOUSEKEEPING:

@testset "local_structured_3d_grids.jl" begin



    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for tt in types

        # setup types and inputs
        t = tt <: Integer ? Float64 : tt <: BigInt ? BigFloat : tt
        N_x, N_y, N_z = rand(2:50, 3)

        #=====================================================================#
        #::. LocalStructured3DGrid
        #=====================================================================#
        x1, y1, z1 = rand(-100:0, 3) |> Vector{t}
        x2, y2, z2 = rand(1:100, 3)  |> Vector{t}
        grid = LocalStructured3DGrid((x1, x2), (y1, y2), (z1, z2), N_x, N_y, N_z)

        @test typeof(grid) <: LocalStructured3DGrid{t}
        @test length(grid) == N_x*N_y*N_z == length(coords(grid))
        @test size(grid) == (N_x, N_y, N_z)

        # surface coords test
        scoords = surfacecoords(grid)
        for sc in scoords; @test isapprox(sc.z, z1 + (z2-z1)/N_z/2); end

        # areas and volume tests
        if t != Float16 # Float16's domain is too small >> Inf16...
            A = (x2-x1) * (y2-y1)
            @test isapprox(sum(areas(grid))/grid.N_z, A)

            V = (x2-x1) * (y2-y1) * (z2-z1)
            @test isapprox(sum(volumes(grid)), V)
        end

        # coord2idx
        for _ in 1:100
            N = rand(1:length(grid))
            coord = coords(grid)[N]
            @test coord2idx(grid, coord) == N
        end
        @test coord2idx(grid, x1-1, y1, z1) == 0


        #=====================================================================#
        #::. LocalStructured3DGrid_Exponential
        #=====================================================================#
        c = rand()*10
        grid = LocalStructured3DGrid_Exponential((x1, x2), (y1, y2), (z1, z2), N_x, N_y, N_z; c=c)

        @test typeof(grid) <: LocalStructured3DGrid_Exponential{t}
        @test length(grid) == N_x*N_y*N_z == length(coords(grid))
        @test size(grid) == (N_x, N_y, N_z)

        # surface coords test
        scoords = surfacecoords(grid)
        zmin = min([vec(coord)[3] for coord in coords(grid)]...)
        for sc in scoords; @test isapprox(sc.z, zmin); end

        # areas and volume tests
        if t != Float16 # Float16's domain is too small >> Inf16...
            A = (x2-x1) * (y2-y1)
            @test isapprox(sum(areas(grid))/grid.N_z, A)

            V = (x2-x1) * (y2-y1) * (z2-z1)
            @test isapprox(sum(volumes(grid)), V)
        end
    end

end

nothing
