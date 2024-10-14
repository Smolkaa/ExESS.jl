print("TESTING: misc > solar_angle.jl")

@testset verbose=true "solar_angle.jl ............." begin

    #::. type-stability
    types = [Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t1 in types, t2 in types

        # setup inputs
        theta = 2*rand()*pi - 1pi
        theta = t1 <: Integer ? round(t1, theta) : t1(theta)
        phi = rand()*pi - pi/2
        phi = t2 <: Integer ? round(t2, phi) : t2(phi)
        xs = GlobalSphericalPosition(Int16(1), theta, phi)
        xc = GlobalCartesianPosition(xs)

        # setup out type
        t_out = promote_type(t1, t2)
        if t_out <: Integer; t_out = promote_type(t_out, Float64); end

        # tests
        @test solar_angle(theta, phi) isa t_out
        @test solar_angle(xs) isa t_out
        @test solar_angle(xc) isa t_out
        @test all(t -> (solar_angle(t, theta, phi) isa t), rand(types[4:7], 1000))
    end


    #::. behaviour
    # if one angle is zero, the other solid angle is the solar angle
    @test all(r -> isapprox(solar_angle(0, r), r; rtol=1e-4), rand(1000))
    @test all(r -> isapprox(solar_angle(r, 0), r; rtol=1e-4), rand(1000))

    # if one angle is pi/2, the solar angle is pi/2
    @test all(r -> isapprox(solar_angle(r, pi/2), pi/2; rtol=1e-4), rand(1000)*pi .- pi/2)
    @test all(r -> isapprox(solar_angle(pi/2, r), pi/2; rtol=1e-4), rand(1000)*2*pi .- pi)

    # testing special positions & types
    @test solar_angle(pi, pi) == 0

end

println("\rTESTING: misc > solar_angle.jl - DONE")
nothing
