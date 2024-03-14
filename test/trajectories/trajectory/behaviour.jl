using LinearAlgebra

@testset "behaviour" begin

    # Monte Carlo test for correct starting and landing behaviour
    for _ in 1:100
        # setup
        x0 = rand(GlobalSphericalPosition, SolarSurfaceDistribution(LUNAR_RADIUS))
        v0 = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(250, amu2kg(4)))
        traj = trajectory(x0, v0, ddx_gravity)

        # test launch
        @test isapprox(ExESS._get(traj(0)[1:3]), ExESS._get(GlobalCartesianVelocity(x0, v0)); rtol=1e-6)
        @test isapprox(ExESS._get(traj(0)[4:6]), ExESS._get(GlobalCartesianPosition(x0)); rtol=1e-6)
        @test isapprox(norm(traj(0)[4:6]), LUNAR_RADIUS; rtol=1e-6)
        
        # test landing/escape
        if norm(v0) < escape_velocity(LUNAR_RADIUS, LUNAR_MASS)
            @test isapprox(norm(traj[end][1:3]), norm(v0); rtol=1e-4)
            @test isapprox(norm(traj[end][4:6]), LUNAR_RADIUS; rtol=1e-6)
        else
            @test isapprox(norm(traj[end][4:6]), 1e9, rtol=1e-6)
        end
    end
    
end

nothing