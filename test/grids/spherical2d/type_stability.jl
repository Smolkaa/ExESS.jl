@testset "type stability" begin

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
    end
end