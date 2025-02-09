@testset "type stability - subtypes" begin

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

    # velocity vectors
    for t1 in types, t2 in types
        t_out_1 = t1 <: Integer ? promote_type(t1, Float64) : t1
        t_out_2 = t2 <: Integer ? promote_type(t2, Float64) : t2
        t_out_12 = promote_type(t_out_1, t_out_2)

        x = Tuple(t1 <: Integer ? rand(1:100, 3) |> Vector{t1} : rand(t1, 3) * 100)
        v = Tuple(t2 <: Integer ? rand(1:100, 3) |> Vector{t2} : rand(t2, 3) * 100)
        xc, xlc, xs = GlobalCartesianPosition(x), LocalCartesianPosition(x), GlobalSphericalPosition(x)
        vc, vlc = GlobalCartesianVelocity(v), LocalCartesianVelocity(v)
        @test typeof(vc) == GlobalCartesianVelocity{t_out_2}
        @test typeof(vlc) == LocalCartesianVelocity{t_out_2}

        @test typeof(GlobalCartesianVelocity(vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xc, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xlc, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xs, vc)) == GlobalCartesianVelocity{t_out_2}
        @test typeof(GlobalCartesianVelocity(xc, vlc)) == GlobalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(GlobalCartesianVelocity(xlc, vlc))
        @test typeof(GlobalCartesianVelocity(xs, vlc)) == GlobalCartesianVelocity{t_out_12}

        @test typeof(LocalCartesianVelocity(xc, vc)) == LocalCartesianVelocity{t_out_12}
        @test_throws MethodError typeof(LocalCartesianVelocity(xlc, vc))
        @test typeof(LocalCartesianVelocity(xs, vc)) == LocalCartesianVelocity{t_out_12}
        @test typeof(LocalCartesianVelocity(vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xc, vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xlc, vlc)) == LocalCartesianVelocity{t_out_2}
        @test typeof(LocalCartesianVelocity(xs, vlc)) == LocalCartesianVelocity{t_out_2}

    end

end
