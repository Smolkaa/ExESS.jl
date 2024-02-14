@testset "base extensions - velocity" begin

    a, v, vr = rand(), rand(3), rand(3) .* [1, 2pi, pi] .- [0, pi, pi/2]
    vt, vrt = Tuple(v), Tuple(vr)
    vc1, vc2 = GlobalCartesianVelocity(rand(3)), GlobalCartesianVelocity(rand(3))
    vs1, vs2 = LocalCartesianVelocity(rand(3)), LocalCartesianVelocity(rand(3))
    vr1, vr2 = GlobalSphericalVelocity(rand(3)), GlobalSphericalVelocity(rand(3))

    @test typeof(vc1 + vt) <: GlobalCartesianVelocity
    @test typeof(vc1 + v) <: GlobalCartesianVelocity
    @test typeof(vt + vc2) <: GlobalCartesianVelocity
    @test typeof(v + vc2) <: GlobalCartesianVelocity
    @test typeof(vc1 + vc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 + vs2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 + vr2) <: GlobalCartesianVelocity

    @test typeof(vs1 + vt) <: LocalCartesianVelocity
    @test typeof(vs1 + v) <: LocalCartesianVelocity
    @test typeof(vt + vs2) <: LocalCartesianVelocity
    @test typeof(v + vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 + vc2) <: LocalCartesianVelocity
    @test typeof(vs1 + vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 + vr2) <: LocalCartesianVelocity

    @test typeof(vr1 + vrt) <: GlobalSphericalVelocity
    @test typeof(vr1 + vr) <: GlobalSphericalVelocity
    @test typeof(vrt + vr2) <: GlobalSphericalVelocity
    @test typeof(vr + vr2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 + vc2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 + vs2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 + vr2) <: GlobalSphericalVelocity

    @test typeof(vc1 - vt) <: GlobalCartesianVelocity
    @test typeof(vc1 - v) <: GlobalCartesianVelocity
    @test typeof(vt - vc2) <: GlobalCartesianVelocity
    @test typeof(v - vc2) <: GlobalCartesianVelocity
    @test typeof(vc1 - vc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 - vs2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 - vr2) <: GlobalCartesianVelocity

    @test typeof(vs1 - vt) <: LocalCartesianVelocity
    @test typeof(vs1 - v) <: LocalCartesianVelocity
    @test typeof(vt - vs2) <: LocalCartesianVelocity
    @test typeof(v - vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 - vc2) <: LocalCartesianVelocity
    @test typeof(vs1 - vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 - vr2) <: LocalCartesianVelocity

    @test typeof(vr1 - vrt) <: GlobalSphericalVelocity
    @test typeof(vr1 - vr) <: GlobalSphericalVelocity
    @test typeof(vrt - vr2) <: GlobalSphericalVelocity
    @test typeof(vr - vr2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 - vc2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 - vs2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 - vr2) <: GlobalSphericalVelocity

    @test typeof(vc1 * a) <: GlobalCartesianVelocity
    @test typeof(a * vc1) <: GlobalCartesianVelocity
    @test typeof(vc1 * vt) <: GlobalCartesianVelocity
    @test typeof(vc1 * v) <: GlobalCartesianVelocity
    @test typeof(vt * vc2) <: GlobalCartesianVelocity
    @test typeof(v * vc2) <: GlobalCartesianVelocity
    @test typeof(vc1 * vc2) <: Real
    @test_throws MethodError typeof(vc1 * vs2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 * vr2) <: GlobalCartesianVelocity

    @test typeof(vs1 * a) <: LocalCartesianVelocity
    @test typeof(a * vs1) <: LocalCartesianVelocity
    @test typeof(vs1 * vt) <: LocalCartesianVelocity
    @test typeof(vs1 * v) <: LocalCartesianVelocity
    @test typeof(vt * vs2) <: LocalCartesianVelocity
    @test typeof(v * vs2) <: LocalCartesianVelocity
    @test typeof(vs1 * vs2) <: Real
    @test_throws MethodError typeof(vs1 * vc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 * vr2) <: LocalCartesianVelocity

    @test typeof(vr1 * a) <: GlobalSphericalVelocity
    @test typeof(a * vr1) <: GlobalSphericalVelocity
    @test typeof(vr1 * vrt) <: GlobalSphericalVelocity
    @test typeof(vr1 * vr) <: GlobalSphericalVelocity
    @test typeof(vrt * vr2) <: GlobalSphericalVelocity
    @test typeof(vr * vr2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 * vc2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 * vs2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 * vr2) <: Real

    @test typeof(vc1 / a) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(a / vc1) <: GlobalCartesianVelocity
    @test typeof(vc1 / vt) <: GlobalCartesianVelocity
    @test typeof(vc1 / v) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vt / vc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(v / vc2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 / vc2) <: Real
    @test_throws MethodError typeof(vc1 / vs2) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(vc1 / vr2) <: GlobalCartesianVelocity

    @test typeof(vs1 / a) <: LocalCartesianVelocity
    @test_throws MethodError typeof(a / vs1) <: LocalCartesianVelocity
    @test typeof(vs1 / vt) <: LocalCartesianVelocity
    @test typeof(vs1 / v) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vt / vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(v / vs2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 / vc2) <: LocalCartesianVelocity
    @test_throws MethodError typeof(vs1 / vs2) <: Real
    @test_throws MethodError typeof(vs1 / vr2) <: LocalCartesianVelocity

    @test typeof(vr1 / a) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(a / vr1) <: GlobalSphericalVelocity
    @test typeof(vr1 / vrt) <: GlobalSphericalVelocity
    @test typeof(vr1 / vr) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vrt / vr2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr / vr2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 / vc2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 / vs2) <: GlobalSphericalVelocity
    @test_throws MethodError typeof(vr1 / vr2) <: Real

end