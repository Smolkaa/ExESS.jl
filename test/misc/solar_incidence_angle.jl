print("TESTING: misc > solar_incidence_angle.jl")

@testset verbose=true "solar_incidence_angle.jl ................." begin

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
        @test solar_incidence_angle(theta, phi) isa t_out
        @test solar_incidence_angle(xs) isa t_out
        @test solar_incidence_angle(xc) isa t_out
        @test all(t -> (solar_incidence_angle(t, theta, phi) isa t), rand(types[4:7], 1000))
    end


    #::. behaviour
    # if one angle is zero, the other solid angle is the solar angle
    @test all(r -> isapprox(solar_incidence_angle(0, r), r; rtol=1e-4), rand(1000))
    @test all(r -> isapprox(solar_incidence_angle(r, 0), r; rtol=1e-4), rand(1000))

    # if one angle is pi/2, the solar angle is pi/2
    @test all(r -> isapprox(solar_incidence_angle(r, pi/2), pi/2; rtol=1e-4), rand(1000)*pi .- pi/2)
    @test all(r -> isapprox(solar_incidence_angle(pi/2, r), pi/2; rtol=1e-4), rand(1000)*2*pi .- pi)

    # testing special positions & types
    @test solar_incidence_angle(pi, pi) == 0

    # subsolar point: local angle equals slope
    @test all(s -> isapprox(local_solar_incidence_angle(0,0,s,rand()), s; rtol=1e-4), rand(1000)*pi)

    # correct rotation effects
    @test all(t -> isapprox(local_solar_incidence_angle(-t, 0, t, 0), 0; atol=1e-4), rand(1000)*2pi .- pi)
    @test all(p -> isapprox(local_solar_incidence_angle(0, p, p, -pi/2), 0; atol=1e-4), rand(1000)*pi/2)
    @test all(p -> isapprox(local_solar_incidence_angle(0, -p, p, pi/2), 0; atol=1e-4), rand(1000)*pi/2)

    # no slope: local angle equals global angle
    @test all(_ -> begin
        t, p, r = (rand(3) * 2pi .- pi) .* [1, 0.5, 0.5]
        return isapprox(solar_incidence_angle(t, p), local_solar_incidence_angle(t, p, 0, r); rtol=1e-4)
    end, 1:1000)

    # TODO: decl?

end

println("\rTESTING: misc > solar_incidence_angle.jl - DONE")
nothing
