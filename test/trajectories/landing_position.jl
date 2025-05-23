using LinearAlgebra

@testset verbose=true "landing_position.jl ......................" begin

    print("TESTING: trajectories > landing_position.jl")

    #=====================================================================#
    #::. type stability
    #=====================================================================#
    TYPES = [Int32, Int64, BigInt, Float32, Float64, BigFloat]
    for t1 in TYPES, t2 in TYPES, _ in 1:25

        # setup inputs
        x0 = rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS))
        v0 = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
        x0_v = t1 <: Integer ? ceil.(vec(x0)) |> Vector{t1} : t1.(vec(x0))
        v0_v = t2 <: Integer ? ceil.(vec(v0)) |> Vector{t2} : t2.(vec(v0))
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
        if !isnan(Tuple(landing_position(x0_t, v0_t))[1])
            @test landing_position(x0_t, v0_t) isa Tuple{t_out, t_out, t_out}
            @test time_of_flight(x0_t, v0_t) isa t_out
        end
        if !isnan(Tuple(landing_position(x0_v, v0_v))[1])
            @test landing_position(x0_v, v0_v) isa Tuple{t_out, t_out, t_out}
            @test time_of_flight(x0_v, v0_v) isa t_out
        end
        if !isnan(Tuple(landing_position(x0_gc, v0_lc))[1])
            @test landing_position(x0_gc, v0_lc) isa GlobalSphericalPosition{t_out_int}
            @test time_of_flight(x0_gc, v0_lc) isa t_out_int
        end
        if !isnan(Tuple(landing_position(x0_gs, v0_lc))[1])
            @test landing_position(x0_gs, v0_lc) isa GlobalSphericalPosition{t_out_int}
            @test time_of_flight(x0_gs, v0_lc) isa t_out_int
        end
    end



    #=====================================================================#
    #::. behaviour
    #=====================================================================#
    # MC tests for accuracy / correctness
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

    # zero-velocity should land at launch position without throwing NaNs @ zero time-of-flight
    @test all(x0 -> isapprox(landing_position(x0, LocalCartesianVelocity(0, 0, 0)), x0; rtol=1e-4), rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS), 1000))
    @test all(x0 -> isapprox(time_of_flight(x0, LocalCartesianVelocity(0, 0, 0)), 0; rtol=1e-4), rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS), 1000))

    # vertically upward velocity should land at launch position without throwing NaNs
    @test all(x0 -> isapprox(landing_position(x0, LocalCartesianVelocity(0, 0, rand())), x0; rtol=1e-4), rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS), 1000))

    # horizontal velocity should land at launch position without throwing NaNs @ zero time-of-flight
    @test all(x0 -> isapprox(landing_position(x0, LocalCartesianVelocity(rand(), rand(), 0)), x0; rtol=1e-4), rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS), 1000))
    @test all(x0 -> isapprox(time_of_flight(x0, LocalCartesianVelocity(rand(), rand(), 0)), 0; rtol=1e-4), rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS), 1000))




    println("\rTESTING: trajectories > landing_position.jl - DONE")
end

nothing
