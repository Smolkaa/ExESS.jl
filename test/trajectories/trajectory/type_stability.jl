#=
    Note that 16-bit types are not tested here as they are usually not accurate enough for
    typical trajectories. Interpolation based on 16-bit types is not supported and, thus,
    will be promoted.
=#
@testset "type stability" begin

    types = [Int32, Int64, BigInt, Float32, Float64, BigFloat]
    for t1 in types, t2 in types

        # setup inputs - position
        x0_t = rand(SolarSurfaceDistribution(LUNAR_RADIUS))
        x0_t = t1 <: Integer ? t1.(round.(t1, x0_t)) : t1.(x0_t)
        x0_v = [x0_t...]
        x0_gsp = GlobalSphericalPosition(x0_t)
        x0_gcp = GlobalCartesianPosition(x0_gsp)

        # setup inputs - velocity
        v0_t = rand(MBFluxVelocityDistribution(rand()*200, amu2kg(rand()*100)))
        v0_t = t2 <: Integer ? t2.(round.(t2, v0_t)) : t2.(v0_t)
        v0_lcv = LocalCartesianVelocity(v0_t)
        v0_gcv = GlobalCartesianVelocity{eltype(v0_lcv)}(ExESS._get(GlobalCartesianVelocity(x0_gsp, v0_lcv))...)
        v0_v = [v0_t...]

        # prepare output type
        t_out = promote_type(t1, t2)
        if t_out <: Integer; t_out = promote_type(t_out, Float64); end

        t_rand = rand([Float16, Float32, Float64, BigFloat])
        t_out_rand = promote_type(t_out, t_rand)


        # trajectory w/ tuples & vectors
        traj = trajectory(x0_t, v0_t, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalSphericalPosition(traj, rand(t_rand)) isa GlobalSphericalPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}

        traj = trajectory(x0_v, v0_v, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}

        
        # trajectory w/ custom types
        t_out = t1 <: Integer || t2 <: Integer ? promote_type(t_out, Float64) : t_out
        t_out_rand = promote_type(t_out, t_rand)

        traj = trajectory(x0_gsp, v0_lcv, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}

        traj = trajectory(x0_gsp, v0_gcv, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}

        traj = trajectory(x0_gcp, v0_lcv, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}

        traj = trajectory(x0_gcp, v0_gcv, ddx_gravity)
        @test traj[rand(1:length(traj))][1:6] isa Vector{t_out}
        @test traj(rand(t_rand))[1:6] isa Vector{t_out_rand}
        @test GlobalCartesianPosition(traj, rand(t_rand)) isa GlobalCartesianPosition{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test GlobalCartesianVelocity(traj, rand(t_rand)) isa GlobalCartesianVelocity{t_out_rand}
        @test LocalCartesianVelocity(traj, rand(t_rand)) isa LocalCartesianVelocity{t_out_rand}
    end
end

nothing