@testset verbose=true "lunar_surface_temperatures_CRIDER2002.jl" begin

    ########################################################################################
    #::. TYPE STABILITY
    ########################################################################################
    TYPES = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t1 in TYPES, t3 in TYPES[5:end]

        # test single input types (sza)
        if t1 in TYPES[1:3]
            @test all(x -> lunar_surface_temperatures_CRIDER2002(x) isa Float64, t1.(rand(-1:1, 1000)))
            @test all(x -> lunar_surface_temperatures_CRIDER2002(t3, x) isa t3, t1.(rand(-1:1, 1000)))
            @test lunar_surface_temperatures_CRIDER2002(t1.(rand(-1:1, 1000))) isa Vector{Float64}
            @test lunar_surface_temperatures_CRIDER2002(t3, t1.(rand(-1:1, 1000))) isa Vector{t3}
        elseif t1 == TYPES[4]
            @test all(x -> lunar_surface_temperatures_CRIDER2002(x) isa BigFloat, t1.(rand(-1:1, 1000)))
            @test all(x -> lunar_surface_temperatures_CRIDER2002(t3, x) isa t3, t1.(rand(-1:1, 1000)))
            @test lunar_surface_temperatures_CRIDER2002(t1.(rand(-1:1, 1000))) isa Vector{BigFloat}
            @test lunar_surface_temperatures_CRIDER2002(t3, t1.(rand(-1:1, 1000))) isa Vector{t3}
        else
            @test all(x -> lunar_surface_temperatures_CRIDER2002(x) isa t1, rand(t1, 1000))
            @test all(x -> lunar_surface_temperatures_CRIDER2002(t3, x) isa t3, rand(t1, 1000))
            @test lunar_surface_temperatures_CRIDER2002(rand(t1, 1000)) isa Vector{t1}
            @test lunar_surface_temperatures_CRIDER2002(t3, rand(t1, 1000)) isa Vector{t3}
        end

        # test two input types (lon, lat)
        for t2 in TYPES
            tout = float(promote_type(t1, t2))
            in1 = t1 <: Integer ? t1(rand(-1:1)) : rand(t1) .* pi/2
            in2 = t2 <: Integer ? t2(rand(-1:1)) : rand(t2) .* pi/2

            @test lunar_surface_temperatures_CRIDER2002(in1, in2) isa tout
            @test lunar_surface_temperatures_CRIDER2002(t3, in1, in2) isa t3
            @test lunar_surface_temperatures_CRIDER2002(repeat([in1],100), repeat([in2],100)) isa Vector{tout}
            @test lunar_surface_temperatures_CRIDER2002(t3, repeat([in1],100), repeat([in2],100)) isa Vector{t3}
        end
    end

    # test grid inputs
    for t in TYPES[5:end]
        grid = Spherical2DGrid(1, rand(10:100), rand(10:50))
        @test lunar_surface_temperatures_CRIDER2002(grid) isa Vector{Float64}
        @test lunar_surface_temperatures_CRIDER2002(t, grid) isa Vector{t}

        grid = Spherical3DGrid([1,2,3,4,5], rand(10:100), rand(10:50))
        @test lunar_surface_temperatures_CRIDER2002(grid) isa Vector{Float64}
        @test lunar_surface_temperatures_CRIDER2002(t, grid) isa Vector{t}
    end

    ########################################################################################
    #::. BEHAVIOR
    ########################################################################################

    # test supported input types
    @test lunar_surface_temperatures_CRIDER2002(rand()) isa Float64
    @test lunar_surface_temperatures_CRIDER2002(rand(10)) isa Vector{Float64}
    @test lunar_surface_temperatures_CRIDER2002(rand(), rand()) isa Float64
    @test lunar_surface_temperatures_CRIDER2002(rand(100), rand(100)) isa Vector{Float64}
    @test lunar_surface_temperatures_CRIDER2002(GlobalCartesianPosition(rand(3))) isa Float64
    @test lunar_surface_temperatures_CRIDER2002(GlobalSphericalPosition(rand(3))) isa Float64
    @test lunar_surface_temperatures_CRIDER2002(Spherical2DGrid(1, rand(10:100), rand(10:50))) isa Vector{Float64}

    # daytime temperature: lower for higher sza (at default kwargs)
    @test all(x -> lunar_surface_temperatures_CRIDER2002(x) > lunar_surface_temperatures_CRIDER2002(1.1*x), rand(1000))

    # nighttime temperature: constant
    @test all(x -> lunar_surface_temperatures_CRIDER2002(x; T0=100) == 100, pi/2 .+ rand(1000))
    @test all(T0 -> lunar_surface_temperatures_CRIDER2002(pi/2+0.1; T0=T0) == T0, rand(1000)*300)

    # same output for sza and (lon, lat) inputs at same angular positions
    @test all(x -> isapprox(lunar_surface_temperatures_CRIDER2002(x), lunar_surface_temperatures_CRIDER2002(0, x); rtol=1e-6), rand(1000))
    @test all(x -> isapprox(lunar_surface_temperatures_CRIDER2002(x), lunar_surface_temperatures_CRIDER2002(x, 0); rtol=1e-6), rand(1000))

    # errors if sza is outside (-π, π)
    @test_throws AssertionError lunar_surface_temperatures_CRIDER2002(pi + rand())
end
