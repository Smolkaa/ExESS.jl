@testset "base extensions - positions" begin
    
    a, x, r = rand(), rand(3), rand(3) .* [1, 2pi, pi] .- [0, pi, pi/2]
    xt, rt = Tuple(x), Tuple(r)
    x1, x2 = GlobalCartesianPosition(rand(3)), GlobalCartesianPosition(rand(3))
    l1, l2 = LocalCartesianPosition(rand(3)), LocalCartesianPosition(rand(3))
    r1, r2 = GlobalSphericalPosition(rand(3)), GlobalSphericalPosition(rand(3))

    @test typeof(x1 + xt) <: GlobalCartesianPosition
    @test typeof(x1 + x) <: GlobalCartesianPosition
    @test typeof(xt + x2) <: GlobalCartesianPosition
    @test typeof(x + x2) <: GlobalCartesianPosition
    @test typeof(x1 + x2) <: GlobalCartesianPosition
    @test_throws MethodError typeof(x1 + l2) <: GlobalCartesianPosition
    @test typeof(x1 + r1) <: GlobalCartesianPosition

    @test typeof(l1 + xt) <: LocalCartesianPosition
    @test typeof(l1 + x) <: LocalCartesianPosition
    @test typeof(xt + l2) <: LocalCartesianPosition
    @test typeof(x + l2) <: LocalCartesianPosition
    @test typeof(l1 + l2) <: LocalCartesianPosition
    @test_throws MethodError typeof(l1 + x2) <: LocalCartesianPosition
    @test_throws MethodError typeof(l1 + r2) <: LocalCartesianPosition

    @test typeof(r1 + rt) <: GlobalSphericalPosition
    @test typeof(r1 + r) <: GlobalSphericalPosition
    @test typeof(rt + r2) <: GlobalSphericalPosition
    @test typeof(r + r2) <: GlobalSphericalPosition
    @test typeof(r1 + x2) <: GlobalSphericalPosition
    @test_throws MethodError typeof(r1 + l2) <: GlobalSphericalPosition
    @test typeof(r1 + r2) <: GlobalSphericalPosition
    
    @test typeof(x1 - xt) <: GlobalCartesianPosition
    @test typeof(x1 - x) <: GlobalCartesianPosition
    @test typeof(xt - x2) <: GlobalCartesianPosition
    @test typeof(x - x2) <: GlobalCartesianPosition
    @test typeof(x1 - x2) <: GlobalCartesianPosition
    @test_throws MethodError typeof(x1 - l2) <: GlobalCartesianPosition
    @test typeof(x1 - r1) <: GlobalCartesianPosition

    @test typeof(l1 - xt) <: LocalCartesianPosition
    @test typeof(l1 - x) <: LocalCartesianPosition
    @test typeof(xt - l2) <: LocalCartesianPosition
    @test typeof(x - l2) <: LocalCartesianPosition
    @test typeof(l1 - l2) <: LocalCartesianPosition
    @test_throws MethodError typeof(l1 - x2) <: LocalCartesianPosition
    @test_throws MethodError typeof(l1 - r2) <: LocalCartesianPosition

    @test typeof(r1 - rt) <: GlobalSphericalPosition
    @test typeof(r1 - r) <: GlobalSphericalPosition
    @test typeof(rt - r2) <: GlobalSphericalPosition
    @test typeof(r - r2) <: GlobalSphericalPosition
    @test typeof(r1 - x2) <: GlobalSphericalPosition
    @test_throws MethodError typeof(r1 - l2) <: GlobalSphericalPosition
    @test typeof(r1 - r2) <: GlobalSphericalPosition

    @test typeof(x1 * a) <: GlobalCartesianPosition
    @test typeof(a * x2) <: GlobalCartesianPosition
    @test typeof(x1 * xt) <: GlobalCartesianPosition
    @test typeof(x1 * x) <: GlobalCartesianPosition
    @test typeof(xt * x2) <: GlobalCartesianPosition
    @test typeof(x * x2) <: GlobalCartesianPosition
    @test typeof(x1 * x2) <: Real
    @test_throws MethodError typeof(x1 * l2) <: Real
    @test typeof(x1 * r2) <: Real

    @test typeof(l1 * a) <: LocalCartesianPosition
    @test typeof(a * l2) <: LocalCartesianPosition
    @test typeof(l1 * xt) <: LocalCartesianPosition
    @test typeof(l1 * x) <: LocalCartesianPosition
    @test typeof(xt * l2) <: LocalCartesianPosition
    @test typeof(x * l2) <: LocalCartesianPosition
    @test typeof(l1 * l2) <: Real
    @test_throws MethodError typeof(l1 * x2) <: Real 
    @test_throws MethodError typeof(l1 * r2) <: Real

    @test typeof(r1 * a) <: GlobalSphericalPosition
    @test typeof(a * r2) <: GlobalSphericalPosition
    @test typeof(r1 * rt) <: GlobalSphericalPosition
    @test typeof(r1 * r) <: GlobalSphericalPosition
    @test typeof(rt * r2) <: GlobalSphericalPosition
    @test typeof(r * r2) <: GlobalSphericalPosition
    @test typeof(r1 * x2) <: Real
    @test_throws MethodError typeof(r1 * l2) <: Real
    @test typeof(r1 * r2) <: Real
    
    @test typeof(x1 / a) <: GlobalCartesianPosition
    @test_throws MethodError typeof(a / x2) <: GlobalCartesianPosition
    @test typeof(x1 / xt) <: GlobalCartesianPosition
    @test typeof(x1 / x) <: GlobalCartesianPosition
    @test_throws MethodError typeof(xt / x2) <: GlobalCartesianPosition
    @test_throws MethodError typeof(x / x2) <: GlobalCartesianPosition
    @test_throws MethodError typeof(x1 / x2) <: Real

    @test typeof(l1 / a) <: LocalCartesianPosition
    @test_throws MethodError typeof(a / l2) <: LocalCartesianPosition
    @test typeof(l1 / xt) <: LocalCartesianPosition
    @test typeof(l1 / x) <: LocalCartesianPosition
    @test_throws MethodError typeof(xt / l2) <: LocalCartesianPosition
    @test_throws MethodError typeof(x / l2) <: LocalCartesianPosition
    @test_throws MethodError typeof(l1 / l2) <: Real

    @test typeof(r1 / a) <: GlobalSphericalPosition
    @test_throws MethodError typeof(a / r2) <: GlobalSphericalPosition
    @test typeof(r1 / rt) <: GlobalSphericalPosition
    @test typeof(r1 / r) <: GlobalSphericalPosition
    @test_throws MethodError typeof(rt / r2) <: GlobalSphericalPosition
    @test_throws MethodError typeof(r / r2) <: GlobalSphericalPosition
    @test_throws MethodError typeof(r1 / r2) <: Real

end