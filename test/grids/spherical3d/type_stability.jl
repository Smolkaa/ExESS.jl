@testset "type stability" begin

    #
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t in types

        # prepare inputs
        r = t <: Integer ? t.(rand(1:100)) : t.(rand() * 100)
        h = t <: Integer ? t.(rand(1:100, rand(2:10))) : t.(rand(rand(2:10)) * 100)
        N_lon, N_lat = rand(2:180), rand(2:90)

        # prepare expected output type
        t_out = t <: Integer ? promote_type(t, Float64) : t

        # create grids
        grid = Spherical3DGrid(r, h, N_lon, N_lat)
        grid_eqsim = Spherical3DGrid_EqSim(r, h, N_lon, N_lat)
        grid_reduced = Spherical3DGrid_Reduced(r, h, N_lat)
        grid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(r, h, N_lat)

        # create partial grids with random ranges
        lonmin, latmin = rand(2) .* rand([-1,1])
        lonmax, latmax = lonmin + rand()*pi/10, latmin + rand()*pi/10
        pgrid = Spherical3DGrid(r, h, N_lon, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
        pgrid_eqsim = Spherical3DGrid_EqSim(r, h, N_lon, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))
        pgrid_reduced = Spherical3DGrid_Reduced(r, h, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
        pgrid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(r, h, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))

        # test grid type
        @test grid isa Spherical3DGrid{t_out}
        @test grid_eqsim isa Spherical3DGrid_EqSim{t_out}
        @test grid_reduced isa Spherical3DGrid_Reduced{t_out}
        @test grid_reduced_eqsim isa Spherical3DGrid_Reduced_EqSim{t_out}

        # test partial grid type
        @test pgrid isa Spherical3DGrid{t_out}
        @test pgrid_eqsim isa Spherical3DGrid_EqSim{t_out}
        @test pgrid_reduced isa Spherical3DGrid_Reduced{t_out}
        @test pgrid_reduced_eqsim isa Spherical3DGrid_Reduced_EqSim{t_out}

        # test type override
        let tt = rand([Float16, Float32, Float64, BigFloat])
            tgrid = Spherical3DGrid(tt, r, h, N_lon, N_lat)
            tgrid_eqsim = Spherical3DGrid_EqSim(tt, r, h, N_lon, N_lat)
            tgrid_reduced = Spherical3DGrid_Reduced(tt, r, h, N_lat)
            tgrid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(tt, r, h, N_lat)

            tpgrid = Spherical3DGrid(tt, r,h,  N_lon, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
            tpgrid_eqsim = Spherical3DGrid_EqSim(tt, r, h, N_lon, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))
            tpgrid_reduced = Spherical3DGrid_Reduced(tt, r, h, N_lat; lonrange=(lonmin, lonmax), latrange=(latmin, latmax))
            tpgrid_reduced_eqsim = Spherical3DGrid_Reduced_EqSim(tt, r, h, N_lat; lonrange=(lonmin, lonmax), latmax=abs(latmax))

            @test tgrid isa Spherical3DGrid{tt}
            @test tgrid_eqsim isa Spherical3DGrid_EqSim{tt}
            @test tgrid_reduced isa Spherical3DGrid_Reduced{tt}
            @test tgrid_reduced_eqsim isa Spherical3DGrid_Reduced_EqSim{tt}

            @test tpgrid isa Spherical3DGrid{tt}
            @test tpgrid_eqsim isa Spherical3DGrid_EqSim{tt}
            @test tpgrid_reduced isa Spherical3DGrid_Reduced{tt}
            @test tpgrid_reduced_eqsim isa Spherical3DGrid_Reduced_EqSim{tt}
        end
    end
end
