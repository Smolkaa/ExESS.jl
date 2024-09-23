@testset "utility functions" begin

    # azimuth
    @test all(t -> azimuth(t, rand(3)) isa t, [Float16, Float32, Float64, BigFloat])
    for x in [(1,0,0), (1,1,pi/4), (0,1,pi/2), (-1,1,3pi/4), (-1,0,pi), (-1,-1,-3pi/4), (0,-1,-pi/2), (1,-1,-pi/4)]
        @test all(r -> isapprox(azimuth((x[1],x[2],r)), x[3]; rtol=1e-10), rand(100))
        @test all(r -> isapprox(azimuth([x[1],x[2],r]), x[3]; rtol=1e-10), rand(100))
        @test all(r -> isapprox(azimuth(GlobalCartesianPosition(x[1],x[2],r)), x[3]; rtol=1e-10), rand(100))
        @test all(r -> isapprox(azimuth(GlobalCartesianVelocity(x[1],x[2],r)), x[3]; rtol=1e-10), rand(100))
    end

    # elevation

    # rotate_x
    @test rotate_x((1, 0, 0), rand()) isa Tuple
    @test rotate_x([1, 0, 0], rand()) isa Vector
    @test rotate_x(GlobalCartesianPosition(1, 0, 0), rand()) isa GlobalCartesianPosition
    @test rotate_x(LocalCartesianPosition(1, 0, 0), rand()) isa LocalCartesianPosition
    @test rotate_x(GlobalCartesianVelocity(1, 0, 0), rand()) isa GlobalCartesianVelocity
    @test rotate_x(LocalCartesianVelocity(1, 0, 0), rand()) isa LocalCartesianVelocity

    @test all(r -> (rotate_x((1,0,0), r) == (1,0,0)), rand(1000))
    @test isapprox(rotate_x((0,1,0),  pi/2),  (0,0,1); atol=1e-10)
    @test isapprox(rotate_x((0,1,0),    pi), (0,-1,0); atol=1e-10)
    @test isapprox(rotate_x((0,1,0), -pi/2), (0,0,-1); atol=1e-10)
    @test isapprox(rotate_x((0,0,1),  pi/2), (0,-1,0); atol=1e-10)
    @test isapprox(rotate_x((0,0,1),    pi), (0,0,-1); atol=1e-10)
    @test isapprox(rotate_x((0,0,1), -pi/2),  (0,1,0); atol=1e-10)

    # rotat_y
    @test rotate_y((1, 0, 0), rand()) isa Tuple
    @test rotate_y([1, 0, 0], rand()) isa Vector
    @test rotate_y(GlobalCartesianPosition(1, 0, 0), rand()) isa GlobalCartesianPosition
    @test rotate_y(LocalCartesianPosition(1, 0, 0), rand()) isa LocalCartesianPosition
    @test rotate_y(GlobalCartesianVelocity(1, 0, 0), rand()) isa GlobalCartesianVelocity
    @test rotate_y(LocalCartesianVelocity(1, 0, 0), rand()) isa LocalCartesianVelocity

    @test isapprox(rotate_y((1,0,0), pi/2), (0,0,-1); atol=1e-10)
    @test isapprox(rotate_y((1,0,0),   pi), (-1,0,0); atol=1e-10)
    @test isapprox(rotate_y((1,0,0),-pi/2),  (0,0,1); atol=1e-10)
    @test all(r -> (rotate_y((0,1,0), r) == (0,1,0)), rand(1000))
    @test isapprox(rotate_y((0,0,1), pi/2),  (1,0,0); atol=1e-10)
    @test isapprox(rotate_y((0,0,1),   pi), (0,0,-1); atol=1e-10)
    @test isapprox(rotate_y((0,0,1),-pi/2), (-1,0,0); atol=1e-10)

    # rotate_z
    @test rotate_z((1, 0, 0), rand()) isa Tuple
    @test rotate_z([1, 0, 0], rand()) isa Vector
    @test rotate_z(GlobalCartesianPosition(1, 0, 0), rand()) isa GlobalCartesianPosition
    @test rotate_z(LocalCartesianPosition(1, 0, 0), rand()) isa LocalCartesianPosition
    @test rotate_z(GlobalCartesianVelocity(1, 0, 0), rand()) isa GlobalCartesianVelocity
    @test rotate_z(LocalCartesianVelocity(1, 0, 0), rand()) isa LocalCartesianVelocity

    @test isapprox(rotate_z((1,0,0), pi/2), (0,1,0); atol=1e-10)
    @test isapprox(rotate_z((1,0,0),   pi), (-1,0,0); atol=1e-10)
    @test isapprox(rotate_z((1,0,0),-pi/2), (0,-1,0); atol=1e-10)
    @test isapprox(rotate_z((0,1,0), pi/2), (-1,0,0); atol=1e-10)
    @test isapprox(rotate_z((0,1,0),   pi), (0,-1,0); atol=1e-10)
    @test isapprox(rotate_z((0,1,0),-pi/2),  (1,0,0); atol=1e-10)
    @test all(r -> (rotate_z((0,0,1), r) == (0,0,1)), rand(1000))

    # speed

    # zenith

end
