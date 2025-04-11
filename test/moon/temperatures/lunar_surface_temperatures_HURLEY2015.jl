@testset verbose=true "lunar_surface_temperatures_HURLEY2015.jl" begin

    ########################################################################################
    #::. TYPE STABILITY
    ########################################################################################
    TYPES = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for t in TYPES[5:end]

        # positional input
        for t1 in TYPES, t2 in TYPES
            tout = float(promote_type(t1, t2))
            in1 = t1 <: Integer ? t1(rand(-1:1)) : rand(t1) .* pi/2
            in2 = t2 <: Integer ? t2(rand(-1:1)) : rand(t2) .* pi/2

            @test lunar_surface_temperatures_HURLEY2015(in1, in2) isa tout
            @test lunar_surface_temperatures_HURLEY2015(t, in1, in2) isa t
            @test lunar_surface_temperatures_HURLEY2015(repeat([in1],100), repeat([in2],100)) isa Vector{tout}
            @test lunar_surface_temperatures_HURLEY2015(t, repeat([in1],100), repeat([in2],100)) isa Vector{t}
        end

        # grid input
        grid = Spherical2DGrid(1, rand(10:100), rand(10:50))
        @test lunar_surface_temperatures_HURLEY2015(grid) isa Vector{Float64}
        @test lunar_surface_temperatures_HURLEY2015(t, grid) isa Vector{t}

        grid = Spherical3DGrid([1,2,3,4,5], rand(10:100), rand(10:50))
        @test lunar_surface_temperatures_HURLEY2015(grid) isa Vector{Float64}
        @test lunar_surface_temperatures_HURLEY2015(t, grid) isa Vector{t}
    end

    ########################################################################################
    #::. BEHAVIOR
    ########################################################################################

    # test all supported input types
    @test lunar_surface_temperatures_HURLEY2015(rand(), rand()) isa Float64
    @test lunar_surface_temperatures_HURLEY2015(rand(10), rand(10)) isa Vector{Float64}
    @test lunar_surface_temperatures_HURLEY2015(GlobalCartesianPosition(rand(3))) isa Float64
    @test lunar_surface_temperatures_HURLEY2015(GlobalSphericalPosition(rand(3))) isa Float64
    @test lunar_surface_temperatures_HURLEY2015(Spherical2DGrid(1, rand(10:100), rand(10:50))) isa Vector{Float64}

    # errors if latitude is out of bounds
    @test_throws AssertionError lunar_surface_temperatures_HURLEY2015(rand(), pi/2+rand())
end
