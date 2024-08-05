using LinearAlgebra


@testset verbose=true "landing_position.jl ........" begin

    print("TESTING: trajectories > landing_position.jl")

    #=====================================================================#
    #::. type stability
    #=====================================================================#
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t1 in types, t2 in types, _ in 1:10

        # setup inputs
        x0_v = t1 <: Integer ? rand(10_000:20_000, 3) |> Vector{t1} : 10_000 .+ rand(t1, 3) * 10_000
        v0_v = t2 <: Integer ? rand(1:20, 3) |> Vector{t2} : rand(t2, 3) * 20
        x0_t = Tuple(x0_v)
        v0_t = Tuple(v0_v)
        x0_gc = GlobalCartesianPosition(x0_t)
        x0_gs = GlobalSphericalPosition(x0_t)
        v0_lc = LocalCartesianVelocity(v0_t)


        # prepare output type
        t_out = promote_type(t1, t2)
        if t_out <: Integer; t_out = promote_type(t_out, Float64); end
        t_out_int = t1 <: Integer || t2 <: Integer ? promote_type(t_out, Float64) : t_out


        # check main functions
        if !isnan(ExESS._get(landing_position(x0_t, v0_t))[1])
            @test typeof(landing_position(x0_t, v0_t)) == Tuple{t_out, t_out, t_out}
            @test typeof(time_of_flight(x0_t, v0_t)) == t_out
        end
        if !isnan(ExESS._get(landing_position(x0_v, v0_v))[1])
            @test typeof(landing_position(x0_v, v0_v)) == Tuple{t_out, t_out, t_out}
            @test typeof(time_of_flight(x0_v, v0_v)) == t_out
        end
        if !isnan(ExESS._get(landing_position(x0_gc, v0_lc))[1])
            @test typeof(landing_position(x0_gc, v0_lc)) == GlobalSphericalPosition{t_out_int}
            @test typeof(time_of_flight(x0_gc, v0_lc)) == t_out_int
        end
        if !isnan(ExESS._get(landing_position(x0_gs, v0_lc))[1])
            @test typeof(landing_position(x0_gs, v0_lc)) == GlobalSphericalPosition{t_out_int}
            @test typeof(time_of_flight(x0_gs, v0_lc)) == t_out_int
        end
    end



    #=====================================================================#
    #::. behaviour (MC test)
    #=====================================================================#
    for _ in 1:2_500
        # inputs
        x0 = rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS))
        v0 = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))

        # solve 3D ODE trajectory
        traj = trajectory(x0, v0, ddx_gravity; reltol=1e-12)

        # calulate landing positions
        land_pos = landing_position(x0, v0)
        land_pos2 = GlobalSphericalPosition(GlobalCartesianPosition(traj.u[end][4:6]))

        # calculate time of flight
        tof = time_of_flight(x0, v0)
        tof2 = traj.t[end]

        # compare landing positions and time of flight if landing occured
        if !isnan(land_pos.r)
            # second condition is to catch bad ODE solver cases too close to v_esc
            if escape_velocity(LUNAR_RADIUS, LUNAR_MASS) - norm(v0) > 10
                @test isapprox(land_pos, land_pos2, rtol=1e-5)

                # note: sometimes the relative tolerance cannot be met, if the tof is too short
                @test isapprox(tof, tof2, rtol=1e-2)
            end
        else
            @test norm(land_pos2) > LUNAR_RADIUS
        end
    end


    println("\rTESTING: trajectories > landing_position.jl - DONE")
end

nothing
