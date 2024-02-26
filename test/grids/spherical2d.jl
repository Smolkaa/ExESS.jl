@testset "spherical2d.jl" begin

    # Monte Carlo tests (N=100)
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:20, t in types

        # prepare inputs
        r = t <: Integer ? t.(rand(1:100)) : t.(rand() * 100)
        N_theta, N_phi = rand(2:180), rand(2:90)

        # prepare expected output type
        t_out = t <: Integer ? promote_type(t, Float64) : t

        # create grids
        grid = Spherical2DGrid(r, N_theta, N_phi)
        grid_eqsim = Spherical2DGrid_EqSim(r, N_theta, N_phi)
        grid_reduced = Spherical2DGrid_Reduced(r, N_phi)
        grid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(r, N_phi)

        # test grid type
        @test grid isa Spherical2DGrid{t_out}
        @test grid_eqsim isa Spherical2DGrid_EqSim{t_out}
        @test grid_reduced isa Spherical2DGrid_Reduced{t_out}
        @test grid_reduced_eqsim isa Spherical2DGrid_Reduced_EqSim{t_out}

        # test type override
        let tt = rand([Float16, Float32, Float64, BigFloat])
            tgrid = Spherical2DGrid(tt, r, N_theta, N_phi)
            tgrid_eqsim = Spherical2DGrid_EqSim(tt, r, N_theta, N_phi)
            tgrid_reduced = Spherical2DGrid_Reduced(tt, r, N_phi)
            tgrid_reduced_eqsim = Spherical2DGrid_Reduced_EqSim(tt, r, N_phi)

            @test tgrid isa Spherical2DGrid{tt}
            @test tgrid_eqsim isa Spherical2DGrid_EqSim{tt}
            @test tgrid_reduced isa Spherical2DGrid_Reduced{tt}
            @test tgrid_reduced_eqsim isa Spherical2DGrid_Reduced_EqSim{tt}
        end

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
        RTOL = t <: Float16 ? 1e-1 : 1e-2 # adjust relative tolerance
        if !isinf(sum(areas(grid)))
            @test isapprox(sum(areas(grid)), 4*pi*r^2, rtol=RTOL)
        end
        if !isinf(sum(areas(grid_eqsim)))
            @test isapprox(sum(areas(grid_eqsim)), 2*pi*r^2, rtol=RTOL)
        end
        if !isinf(sum(areas(grid_reduced)))
            @test isapprox(sum(areas(grid_reduced)), 4*pi*r^2, rtol=RTOL)
        end
        if !isinf(sum(areas(grid_reduced_eqsim)))
            @test isapprox(sum(areas(grid_reduced_eqsim)), 2*pi*r^2, rtol=RTOL)
        end
    end
end