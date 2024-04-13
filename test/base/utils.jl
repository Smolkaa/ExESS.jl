@testset verbose=true "utils.jl ......................" begin

    #::. type stability
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for type in types
        if type == BigInt
            @test amu2kg(type(1))       isa BigFloat
            @test eV2J(type(1))         isa BigFloat
            @test J2eV(type(1))         isa BigFloat
            @test limit_acos(type(1))   isa BigFloat
            @test lng2LT(type(1))       isa BigFloat
            @test LT2lng(type(1))       isa BigFloat
        elseif type <: Integer
            @test amu2kg(type(1))       isa Float64
            @test eV2J(type(1))         isa Float64
            @test J2eV(type(1))         isa Float64
            @test limit_acos(type(1))   isa Float64
            @test lng2LT(type(1))       isa Float64
            @test LT2lng(type(1))       isa Float64
        else
            @test amu2kg(type(1))       isa type
            @test eV2J(type(1))         isa type  # could break for Float16 ?
            @test J2eV(type(1))         isa type  # could break for Float16 ? 
            @test limit_acos(type(1))   isa type
            @test lng2LT(type(1))       isa type
            @test LT2lng(type(1))       isa type
        end
        @test sgn(type(1)) isa type
    end

    #::. eV <> joule conversion
    @test isapprox(eV2J(1), 1.602176634e-19, rtol=1e-6)
    @test all(r -> isapprox(eV2J(J2eV(r)), r, rtol=1e-6), rand(1000) * 10)

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
    @test all(r -> isapprox(LT2lng(lng2LT(r)), r, rtol=1e-6), rand(1000) .* 2pi .- pi)

    #::. sgn
    @test sgn(0) == 1
    @test all(r -> (sgn(r) == 1), rand(1000))
    @test all(r -> (sgn(r) == -1), -rand(1000))
end

nothing