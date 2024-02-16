@testset "temperatures.jl" begin

    #::. type stability
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t1 in types, t2 in types        
        t3 = rand(types) # not inside third loop for performance issues

        # prepare inputs
        in1 = t1 <: Integer ? rand(-1:1) |> t1 : rand(t1)
        in2 = t2 <: Integer ? rand(-1:1) |> t2 : rand(t2)
        in3 = t3 <: Integer ? rand(0:2) |> t3 : rand(t3) * pi * 2

        vin1  = t1 <: Integer ? rand(-1:1, rand(1:9)) |> Vector{t1} : rand(t1, rand(1:9))
        vin2  = t2 <: Integer ? rand(-1:1, rand(1:9)) |> Vector{t2} : rand(t2, rand(1:9))
        vin2b = t2 <: Integer ? rand(-1:1, length(vin1)) |> Vector{t2} : rand(t2, length(vin1))

        # prepare output type
        t_out_1 = t1 <: Integer ? typeof(promote(in1, 1.0)[1]) : t1
        t_out_2 = typeof(promote(in1, in2)[1])
        if t_out_2 <: Integer; t_out_2 = typeof(promote(one(t_out_2), 1.0)[1]); end
        t_out_3 = typeof(promote(in1, in2, in3)[1])
        if t_out_3 <: Integer; t_out_3 = typeof(promote(one(t_out_3), 1.0)[1]); end


        # tests
        @test typeof(lunar_surface_temperatures_BUTLER1997(in1, in2)) == t_out_2
        @test typeof(lunar_surface_temperatures_BUTLER1997(vin1, vin2)) == Matrix{t_out_2}
        @test typeof(lunar_surface_temperatures_BUTLER1997(vin1, vin2b; matrix=false)) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_BUTLER1997(GlobalSphericalPosition(1, in1, in2))) == t_out_2
        @test typeof(lunar_surface_temperatures_BUTLER1997(GlobalSphericalPosition.(1, vin1, vin2b))) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_BUTLER1997(GlobalStructured2DGrid(t_out_1, 1, rand(2:10, 2)...))) == Vector{t_out_1}

        lng, lat, T = lunar_surface_temperatures_diviner(in1 + t1(1))
        @test typeof(lng) == typeof(lat) == typeof(T) == Vector{t_out_1}
        @test typeof(lunar_surface_temperatures_diviner(in3, in1, in2)) == t_out_3
        @test typeof(lunar_surface_temperatures_diviner(in3, vin1, vin2)) == Matrix{t_out_3}
        @test typeof(lunar_surface_temperatures_diviner(in3, vin1, vin2b; matrix=false)) == Vector{t_out_3}
        @test typeof(lunar_surface_temperatures_diviner(in3, GlobalSphericalPosition(1, in1, in2))) == t_out_3
        @test typeof(lunar_surface_temperatures_diviner(in3, GlobalSphericalPosition.(1, vin1, vin2b))) == Vector{t_out_3}
        @test typeof(lunar_surface_temperatures_diviner(in1 + t1(1), GlobalStructured2DGrid(t_out_1, 1, rand(2:10, 2)...))) == Vector{t_out_1}

        @test typeof(lunar_surface_temperatures_diviner_avg(in1, in2)) == t_out_2
        @test typeof(lunar_surface_temperatures_diviner_avg(vin1, vin2)) == Matrix{t_out_2}
        @test typeof(lunar_surface_temperatures_diviner_avg(vin1, vin2b; matrix=false)) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_diviner_avg(GlobalSphericalPosition(1, in1, in2))) == t_out_2
        @test typeof(lunar_surface_temperatures_diviner_avg(GlobalSphericalPosition.(1, vin1, vin2b))) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_diviner_avg(GlobalStructured2DGrid(t_out_1, 1, rand(2:10, 2)...))) == Vector{t_out_1}

        @test typeof(lunar_surface_temperatures_HURLEY2015(in1, in2)) == t_out_2
        @test typeof(lunar_surface_temperatures_HURLEY2015(vin1, vin2)) == Matrix{t_out_2}
        @test typeof(lunar_surface_temperatures_HURLEY2015(vin1, vin2b; matrix=false)) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_HURLEY2015(GlobalSphericalPosition(1, in1, in2))) == t_out_2
        @test typeof(lunar_surface_temperatures_HURLEY2015(GlobalSphericalPosition.(1, vin1, vin2b))) == Vector{t_out_2}
        @test typeof(lunar_surface_temperatures_HURLEY2015(GlobalStructured2DGrid(t_out_1, 1, rand(2:10, 2)...))) == Vector{t_out_1}
    end

    #::. behavior tests
    @test_throws AssertionError lunar_surface_temperatures_BUTLER1997(pi + rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_BUTLER1997(-pi - rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_BUTLER1997((rand()*2-1)*pi, pi/2 + rand())
    @test_throws AssertionError lunar_surface_temperatures_BUTLER1997((rand()*2-1)*pi, - pi/2 - rand())

    @test_throws AssertionError lunar_surface_temperatures_diviner(2pi + rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner(-rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), pi + rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), -pi - rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), (rand()*2-1)*pi, pi/2 + rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), (rand()*2-1)*pi, - pi/2 - rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), [1]*(pi + rand()), [1]*((rand()*2-1)*pi/2))
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), [1]*(-pi - rand()), [1]*((rand()*2-1)*pi/2))
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), [1]*((rand()*2-1)*pi), [1]*(pi/2 + rand()))
    @test_throws AssertionError lunar_surface_temperatures_diviner(rand(), [1]*((rand()*2-1)*pi), [1]*(- pi/2 - rand()))

    @test_throws AssertionError lunar_surface_temperatures_diviner_avg(pi + rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg(-pi - rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg((rand()*2-1)*pi, pi/2 + rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg((rand()*2-1)*pi, - pi/2 - rand())
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg([1]*(pi + rand()), [1]*((rand()*2-1)*pi/2))
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg([1]*(-pi - rand()), [1]*((rand()*2-1)*pi/2))
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg([1]*((rand()*2-1)*pi), [1]*(pi/2 + rand()))
    @test_throws AssertionError lunar_surface_temperatures_diviner_avg([1]*((rand()*2-1)*pi), [1]*(- pi/2 - rand()))

    @test_throws AssertionError lunar_surface_temperatures_HURLEY2015(pi + rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_HURLEY2015(-pi - rand(), (rand()*2-1)*pi/2)
    @test_throws AssertionError lunar_surface_temperatures_HURLEY2015((rand()*2-1)*pi, pi/2 + rand())
    @test_throws AssertionError lunar_surface_temperatures_HURLEY2015((rand()*2-1)*pi, - pi/2 - rand())
end

nothing