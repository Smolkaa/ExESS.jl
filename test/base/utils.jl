@testset verbose=true "utils.jl ................................." begin
    print("TESTING: base > utils.jl")

    #::. type stability
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for type in types
        if type == BigInt
            @test amu2kg(type(1))       isa BigFloat
            @test eV2J(type(1))         isa BigFloat
            @test eV2kJpmol(type(1))    isa BigFloat
            @test J2eV(type(1))         isa BigFloat
            @test kJpmol2eV(type(1))    isa BigFloat
            @test limit_acos(type(1))   isa BigFloat
            @test lon2LT(type(1))       isa BigFloat
            @test LT2lon(type(1))       isa BigFloat
        elseif type <: Integer
            @test amu2kg(type(1))       isa Float64
            @test eV2J(type(1))         isa Float64
            @test eV2kJpmol(type(1))    isa Float64
            @test J2eV(type(1))         isa Float64
            @test kJpmol2eV(type(1))    isa Float64
            @test limit_acos(type(1))   isa Float64
            @test lon2LT(type(1))       isa Float64
            @test LT2lon(type(1))       isa Float64
        else
            @test amu2kg(type(1))       isa type
            @test eV2J(type(1))         isa type  # could break for Float16 ?
            @test eV2kJpmol(type(1))    isa type
            @test J2eV(type(1))         isa type  # could break for Float16 ?
            @test kJpmol2eV(type(1))    isa type
            @test limit_acos(type(1))   isa type
            @test lon2LT(type(1))       isa type
            @test LT2lon(type(1))       isa type
        end
        @test sgn(type(1)) isa type
    end

    #::. eV <> joule conversion
    @test isapprox(eV2J(1), 1.602176634e-19, rtol=1e-6)
    @test all(r -> isapprox(eV2J(J2eV(r)), r, rtol=1e-6), rand(1000) * 10)

    #::. eV <> kJpmol conversion
    @test isapprox(eV2kJpmol(1), 96.4853365, rtol=1e-6)
    @test all(r -> isapprox(eV2kJpmol(kJpmol2eV(r)), r, rtol=1e-6), rand(1000) * 10)

    #::. longitude <> local time conversion
    @test isapprox(lon2LT(-pi),   0, rtol=1e-6)
    @test isapprox(lon2LT(-pi/2), 6, rtol=1e-6)
    @test isapprox(lon2LT(0),    12, rtol=1e-6)
    @test isapprox(lon2LT(pi/2), 18, rtol=1e-6)
    @test isapprox(lon2LT(pi),   24, rtol=1e-6)
    @test isapprox(LT2lon(0),  -pi,   rtol=1e-6)
    @test isapprox(LT2lon(6),  -pi/2, rtol=1e-6)
    @test isapprox(LT2lon(12),  0,    rtol=1e-6)
    @test isapprox(LT2lon(18),  pi/2, rtol=1e-6)
    @test isapprox(LT2lon(24),  pi,   rtol=1e-6)
    @test all(r -> isapprox(LT2lon(lon2LT(r)), r, rtol=1e-6), rand(1000) .* 2pi .- pi)

    #::. pclamp
    function _test_pclamp()
        x, l, u = rand(3)
        u += l
        if l < x < u
            return isapprox(x, pclamp(x, l, u); atol=1e-6)
        elseif x < l
            while x < l; x += u - l; end
            return isapprox(x, pclamp(x, l, u); atol=1e-6)
        else
            while x > u; x -= u - l; end
            return isapprox(x, pclamp(x, l, u); atol=1e-6)
        end
    end
    @test all(_ -> _test_pclamp(), 1:1000) # test proper clamping
    @test all(x -> pclamp(x...) == pclamp(x[1], x[3], x[2]), [rand(3) for _ in 1:1000]) # symmetry

    #::. sgn
    @test sgn(0) == 1
    @test all(r -> (sgn(r) == 1), rand(1000))
    @test all(r -> (sgn(r) == -1), -rand(1000))



    println("\rTESTING: base > utils.jl - DONE")
end

nothing
