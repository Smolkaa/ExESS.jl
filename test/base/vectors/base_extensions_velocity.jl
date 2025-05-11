@testset "base extensions - velocity" begin

    a, v, vr = rand(), rand(3), rand(3) .* [1, 2pi, pi] .- [0, pi, pi/2]
    vt, vrt = Tuple(v), Tuple(vr)
    vgc1, vgc2 = GlobalCartesianVelocity(rand(3)), GlobalCartesianVelocity(rand(3))
    vlc1, vlc2 = LocalCartesianVelocity(rand(3)), LocalCartesianVelocity(rand(3))

    @test typeof(vgc1 + vt) <: GlobalCartesianVelocity
    @test typeof(vgc1 + v) <: GlobalCartesianVelocity
    @test typeof(vt + vgc2) <: GlobalCartesianVelocity
    @test typeof(v + vgc2) <: GlobalCartesianVelocity
    @test typeof(vgc1 + vgc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vgc1 + vlc2) <: GlobalCartesianVelocity

    @test typeof(vlc1 + vt) <: LocalCartesianVelocity
    @test typeof(vlc1 + v) <: LocalCartesianVelocity
    @test typeof(vt + vlc2) <: LocalCartesianVelocity
    @test typeof(v + vlc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vlc1 + vgc2) <: LocalCartesianVelocity
    @test typeof(vlc1 + vlc2) <: LocalCartesianVelocity

    @test typeof(vgc1 - vt) <: GlobalCartesianVelocity
    @test typeof(vgc1 - v) <: GlobalCartesianVelocity
    @test typeof(vt - vgc2) <: GlobalCartesianVelocity
    @test typeof(v - vgc2) <: GlobalCartesianVelocity
    @test typeof(vgc1 - vgc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vgc1 - vlc2) <: GlobalCartesianVelocity

    @test typeof(vlc1 - vt) <: LocalCartesianVelocity
    @test typeof(vlc1 - v) <: LocalCartesianVelocity
    @test typeof(vt - vlc2) <: LocalCartesianVelocity
    @test typeof(v - vlc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vlc1 - vgc2) <: LocalCartesianVelocity
    @test typeof(vlc1 - vlc2) <: LocalCartesianVelocity

    @test typeof(vgc1 * a) <: GlobalCartesianVelocity
    @test typeof(a * vgc1) <: GlobalCartesianVelocity
    @test typeof(vgc1 * vt) <: GlobalCartesianVelocity
    @test typeof(vgc1 * v) <: GlobalCartesianVelocity
    @test typeof(vt * vgc2) <: GlobalCartesianVelocity
    @test typeof(v * vgc2) <: GlobalCartesianVelocity
    @test typeof(vgc1 * vgc2) <: Real
    @test_throws MethodError typeof(vgc1 * vlc2) <: GlobalCartesianVelocity

    @test typeof(vlc1 * a) <: LocalCartesianVelocity
    @test typeof(a * vlc1) <: LocalCartesianVelocity
    @test typeof(vlc1 * vt) <: LocalCartesianVelocity
    @test typeof(vlc1 * v) <: LocalCartesianVelocity
    @test typeof(vt * vlc2) <: LocalCartesianVelocity
    @test typeof(v * vlc2) <: LocalCartesianVelocity
    @test typeof(vlc1 * vlc2) <: Real
    @test_throws MethodError typeof(vlc1 * vgc2) <: LocalCartesianVelocity

    @test typeof(vgc1 / a) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(a / vgc1) <: GlobalCartesianVelocity
    @test typeof(vgc1 / vt) <: GlobalCartesianVelocity
    @test typeof(vgc1 / v) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vt / vgc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(v / vgc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vgc1 / vgc2) <: Real
    @test_throws MethodError typeof(vgc1 / vlc2) <: GlobalCartesianVelocity

    @test typeof(vlc1 / a) <: LocalCartesianVelocity
    @test_throws MethodError typeof(a / vlc1) <: LocalCartesianVelocity
    @test typeof(vlc1 / vt) <: LocalCartesianVelocity
    @test typeof(vlc1 / v) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vt / vlc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(v / vlc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vlc1 / vgc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vlc1 / vlc2) <: Real

end
