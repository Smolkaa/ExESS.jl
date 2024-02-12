@testset "coordinates.jl" begin

    #=====================================================================#
    #::. Constructors; error and type checking
    #=====================================================================#
    xc = GlobalCartesianPosition(rand(3))
    xlc = LocalCartesianPosition(rand(3))
    xs = GlobalSphericalPosition(rand(3))
    vc = GlobalCartesianVelocity(rand(3))
    vlc = LocalCartesianVelocity(rand(3))
    vs = GlobalSphericalVelocity(rand(3))

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
    @test typeof(GlobalCartesianVelocity(xc, vs)) <: GlobalCartesianVelocity
    @test_throws MethodError typeof(GlobalCartesianVelocity(xlc, vs)) <: GlobalCartesianVelocity
    @test typeof(GlobalCartesianVelocity(xs, vs)) <: GlobalCartesianVelocity

    #...


    #=====================================================================#
    #::. error & type checking (internal types)
    #=====================================================================#
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]

    # position vectors
    for t in types
        t_out = t <: Integer ? promote_type(t, Float64) : t

        x = Tuple(t <: Integer ? t.(rand(1:100, 3)) : rand(t, 3) * 100)
        xc, lc, xs = GlobalCartesianPosition(x), LocalCartesianPosition(x), GlobalSphericalPosition(x)
        @test typeof(xc) == GlobalCartesianPosition{t_out}
        @test typeof(lc) == LocalCartesianPosition{t_out}
        @test typeof(xs) == GlobalSphericalPosition{t_out}
        
        @test typeof(GlobalCartesianPosition(xc)) == GlobalCartesianPosition{t_out}
        @test typeof(GlobalCartesianPosition(xs)) == GlobalCartesianPosition{t_out}
        @test_throws MethodError typeof(GlobalCartesianPosition(lc))

        @test_throws MethodError typeof(LocalCartesianPosition(xc))
        @test typeof(LocalCartesianPosition(lc)) == LocalCartesianPosition{t_out}
        @test_throws MethodError typeof(LocalCartesianPosition(xs))

        @test typeof(GlobalSphericalPosition(xc)) == GlobalSphericalPosition{t_out}
        @test typeof(GlobalSphericalPosition(xs)) == GlobalSphericalPosition{t_out}
        @test_throws MethodError typeof(GlobalSphericalPosition(lc))
    end

    for t1 in types, t2 in types
        t_out_1 = t1 <: Integer ? promote_type(t1, Float64) : t1
        t_out_2 = t2 <: Integer ? promote_type(t2, Float64) : t2
        t_out_12 = promote_type(t_out_1, t_out_2)

        x = Tuple(t1 <: Integer ? rand(1:100, 3) |> Vector{t1} : rand(t1, 3) * 100)
        v = Tuple(t2 <: Integer ? rand(1:100, 3) |> Vector{t2} : rand(t2, 3) * 100)
        xc, xlc, xs = GlobalCartesianPosition(x), LocalCartesianPosition(x), GlobalSphericalPosition(x)
        vc, vlc, vs = GlobalCartesianVelocity(v), LocalCartesianVelocity(v), GlobalSphericalVelocity(v)
        @test typeof(vc) == GlobalCartesianVelocity{t_out_2}
        @test typeof(vlc) == LocalCartesianVelocity{t_out_2}
        @test typeof(vs) == GlobalSphericalVelocity{t_out_2}

        @test typeof(GlobalCartesianVelocity(vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xc, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xlc, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xs, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xc, vlc)) == GlobalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(GlobalCartesianVelocity(xlc, vlc))
        @test typeof(GlobalCartesianVelocity(xs, vlc)) == GlobalCartesianVelocity{t_out_12}
        @test typeof(GlobalCartesianVelocity(xc, vs)) == GlobalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(GlobalCartesianVelocity(xlc, vs))
        @test typeof(GlobalCartesianVelocity(xs, vs)) == GlobalCartesianVelocity{t_out_12}

        @test typeof(LocalCartesianVelocity(xc, vc)) == LocalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(LocalCartesianVelocity(xlc, vc))
        @test typeof(LocalCartesianVelocity(xs, vc)) == LocalCartesianVelocity{t_out_12}
        @test typeof(LocalCartesianVelocity(vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xc, vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xlc, vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xs, vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xc, vs)) == LocalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(LocalCartesianVelocity(xlc, vs))
        @test typeof(LocalCartesianVelocity(xs, vs)) == LocalCartesianVelocity{t_out_12}

        @test typeof(GlobalSphericalVelocity(xc, vc)) == GlobalSphericalVelocity{t_out_12}
        @test_throws MethodError typeof(GlobalSphericalVelocity(xlc, vc))
        @test typeof(GlobalSphericalVelocity(xs, vc)) == GlobalSphericalVelocity{t_out_12}
        @test typeof(GlobalSphericalVelocity(xc, vlc)) == GlobalSphericalVelocity{t_out_12}
        @test_throws MethodError typeof(GlobalSphericalVelocity(xlc, vlc))
        @test typeof(GlobalSphericalVelocity(xs, vlc)) == GlobalSphericalVelocity{t_out_12}
        @test typeof(GlobalSphericalVelocity(vs)) == GlobalSphericalVelocity{t_out_2}
        @test typeof(GlobalSphericalVelocity(xc, vs)) == GlobalSphericalVelocity{t_out_2}
        @test typeof(GlobalSphericalVelocity(xlc, vs)) == GlobalSphericalVelocity{t_out_2}
        @test typeof(GlobalSphericalVelocity(xs, vs)) == GlobalSphericalVelocity{t_out_2}
    end

    #=====================================================================#
    #::. Base extensions; position vectors; error and type checking
    #=====================================================================#
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


    #=====================================================================#
    #::. Base extensions; velocity vectors; error and type checking
    #=====================================================================#
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


    #=====================================================================#
    #::. Conversion Tests; mathematical correctness
    #=====================================================================#
    
    
end

nothing