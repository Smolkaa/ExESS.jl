@testset verbose=true "utils.jl ......................" begin

    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]

    #::. amu2kg, limit_acos
    for type in types
        if type == BigInt
            @test typeof(amu2kg(type(1))) == BigFloat
            @test typeof(eV2J(type(1))) == BigFloat
            @test typeof(J2eV(type(1))) == BigFloat
            @test typeof(limit_acos(type(1))) == BigFloat
            @test typeof(lng2LT(type(1))) == BigFloat
            @test typeof(LT2lng(type(1))) == BigFloat
        elseif type <: Integer
            @test typeof(amu2kg(type(1))) == Float64
            @test typeof(eV2J(type(1))) == Float64
            @test typeof(J2eV(type(1))) == Float64
            @test typeof(limit_acos(type(1))) == Float64
            @test typeof(lng2LT(type(1))) == Float64
            @test typeof(LT2lng(type(1))) == Float64
        else
            @test typeof(amu2kg(type(1))) == type
            @test typeof(eV2J(type(1))) == type  # could break for Float16 ?
            @test typeof(J2eV(type(1))) == type  # could break for Float16 ? 
            @test typeof(limit_acos(type(1))) == type
            @test typeof(lng2LT(type(1))) == type
            @test typeof(LT2lng(type(1))) == type
        end
    end

    #::. longitude <> local time conversion
    @test isapprox(lng2LT(-pi),   0, rtol=1e-6)
    @test isapprox(lng2LT(-pi/2), 6, rtol=1e-6)
    @test isapprox(lng2LT(0),    12, rtol=1e-6)
    @test isapprox(lng2LT(pi/2), 18, rtol=1e-6)
    @test isapprox(lng2LT(pi),   24, rtol=1e-6)
    @test isapprox(LT2lng(0),  -pi,   rtol=1e-6)
    @test isapprox(LT2lng(6),  -pi/2, rtol=1e-6)
    @test isapprox(LT2lng(12),  0,    rtol=1e-6)
    @test isapprox(LT2lng(18),  pi/2, rtol=1e-6)
    @test isapprox(LT2lng(24),  pi,   rtol=1e-6)
    for _ in 1:100; r = rand()*2pi-pi; @test isapprox(LT2lng(lng2LT(r)), r, rtol=1e-6); end
end

nothing