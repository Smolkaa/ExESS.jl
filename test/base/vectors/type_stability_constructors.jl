@testset "type stability - constructors" begin

    xc = GlobalCartesianPosition(rand(3))
    xlc = LocalCartesianPosition(rand(3))
    xs = GlobalSphericalPosition(rand(3))
    vc = GlobalCartesianVelocity(rand(3))
    vlc = LocalCartesianVelocity(rand(3))

    @test typeof(GlobalCartesianPosition(rand(3)...)) <: GlobalCartesianPosition
    @test typeof(GlobalCartesianPosition(rand(3))) <: GlobalCartesianPosition
    @test typeof(GlobalCartesianPosition(rand(3) |> Tuple)) <: GlobalCartesianPosition
    @test typeof(GlobalCartesianPosition(xc)) <: GlobalCartesianPosition
    @test_throws MethodError typeof(GlobalCartesianPosition(xlc)) <: GlobalCartesianPosition
    @test typeof(GlobalCartesianPosition(xs)) <: GlobalCartesianPosition

    @test typeof(LocalCartesianPosition(rand(3)...)) <: LocalCartesianPosition
    @test typeof(LocalCartesianPosition(rand(3))) <: LocalCartesianPosition
    @test typeof(LocalCartesianPosition(rand(3) |> Tuple)) <: LocalCartesianPosition
    @test_throws MethodError typeof(LocalCartesianPosition(xc)) <: LocalCartesianPosition
    @test typeof(LocalCartesianPosition(xlc)) <: LocalCartesianPosition
    @test_throws MethodError typeof(LocalCartesianPosition(xs)) <: LocalCartesianPosition

    @test typeof(GlobalSphericalPosition(rand(3)...)) <: GlobalSphericalPosition
    @test typeof(GlobalSphericalPosition(rand(3))) <: GlobalSphericalPosition
    @test typeof(GlobalSphericalPosition(rand(3) |> Tuple)) <: GlobalSphericalPosition
    @test typeof(GlobalSphericalPosition(xc)) <: GlobalSphericalPosition
    @test_throws MethodError typeof(GlobalSphericalPosition(xlc)) <: GlobalSphericalPosition
    @test typeof(GlobalSphericalPosition(xs)) <: GlobalSphericalPosition

    @test typeof(GlobalCartesianVelocity(rand(3)...)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(rand(3))) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(rand(3) |> Tuple)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(vc)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xc, vc)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xlc, vc)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xs, vc)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xc, vlc)) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(GlobalCartesianVelocity(xlc, vlc)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xs, vlc)) <: GlobalCartesianVelocity

    @test LocalCartesianVelocity(rand(3)...) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(rand(3)) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(rand(3) |> Tuple) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(vlc) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(xc, vc) isa LocalCartesianVelocity
    @test_throws MethodError LocalCartesianVelocity(xlc, vc) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(xs, vc) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(xc, vlc) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(xlc, vlc) isa LocalCartesianVelocity
    @test LocalCartesianVelocity(xs, vlc) isa LocalCartesianVelocity

end
