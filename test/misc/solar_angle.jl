print("TESTING: misc > solar_angle.jl")

@testset verbose=true "solar_angle.jl ............." begin

    #::. type_stability
    types = [Int32, Int64, BigInt, Float32, Float64, BigFloat]
    for t1 in types, t2 in types

        theta = t1 <: Integer ? t1(rand(-3:3)) : t1(rand()*2pi - pi)
        phi = t2 <: Integer ? t2(rand(-1:1)) : t2(rand()*pi - pi/2)
        xs = GlobalSphericalPosition(1, theta, phi)
        x = GlobalCartesianPosition(xs)

        t_out = promote_type(t1, t2)
        if t_out <: Integer; t_out = promote_type(t_out, Float64); end

        @test solar_angle(theta, phi) isa t_out
        @test solar_angle(xs) isa t_out
        @test solar_angle(x) isa t_out
    end

    #::. physical behaviour
    @test all(p -> isapprox(solar_angle(0, p), abs(p); atol=1e-6), rand(1000).*pi .- pi/2)
    @test all(t -> isapprox(solar_angle(t, 0), abs(t); atol=1e-6), rand(1000).*2pi .- pi)
end

println("\rTESTING: misc > solar_angle.jl - DONE")
nothing
