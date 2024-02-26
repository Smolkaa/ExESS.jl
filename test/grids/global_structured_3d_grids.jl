@testset "global_structured_3d_grids.jl" begin

    ########################################################################################
    # GlobalStructured3DGrid                                                               #
    ########################################################################################

    #::. type stability and constructor tests
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:100

        # input types
        t1, t2, t3 = rand(types, 3)
        r0 = t2(rand(100:1000))
        h = t3.(rand(100:1000, rand(2:10)))
        r = r0 .+ accumulate(+,sort(h)) |> Vector{t2}

        # outputs
        if t1 <: Integer
            @test_throws MethodError GlobalStructured3DGrid(t1, r0, h, rand(2:10,2)...)
            @test_throws MethodError GlobalStructured3DGrid(t1, r, rand(2:10,2)...)
        else
            @test GlobalStructured3DGrid(t1, r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid{t1}
            @test GlobalStructured3DGrid(t1, r, rand(2:10,2)...) isa GlobalStructured3DGrid{t1}
        end

        if t2 <: Integer
            if t2 == BigInt
                @test GlobalStructured3DGrid(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid{BigFloat}
                @test GlobalStructured3DGrid(r, rand(2:10,2)...) isa GlobalStructured3DGrid{BigFloat}
            else
                @test GlobalStructured3DGrid(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid{Float64}
                @test GlobalStructured3DGrid(r, rand(2:10,2)...) isa GlobalStructured3DGrid{Float64}
            end
        else
            @test GlobalStructured3DGrid(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid{t2}
            @test GlobalStructured3DGrid(r, rand(2:10,2)...) isa GlobalStructured3DGrid{t2}
        end
    end

    #::. area and volume tests
    r, h = Float64.(rand(2) .* 1000)
    grid = GlobalStructured3DGrid(r, [h], rand(10:20,2)...)
    Nsurface = length(surfacecoords(grid))
    
    @test isapprox(sum(areas(grid)[1:Nsurface]), 4 * pi * r^2; rtol=1e-4)
    @test isapprox(sum(volumes(grid)), 4/3 * pi * ((r+h)^3 - r^3); rtol=1e-4)

    #::. coord2idx
    r, h = rand()*100, sort(rand(10)*100)
    grid = GlobalStructured3DGrid(r, h, rand(10:20,2)...)
    Xs = coords(grid)
    for _ in 1:100
        idx = rand(1:length(Xs))
        @test coord2idx(grid, Xs[idx]) == idx
    end

    #::. surfacecoords
    Xssurf = surfacecoords(grid)
    for xs in Xssurf; @test xs.r == r; end

    
    ########################################################################################
    # GlobalStructured3DGrid_EqSim                                                         #
    ########################################################################################

    #::. type stability and constructor tests
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:100

        # input types
        t1, t2, t3 = rand(types, 3)
        r0 = t2(rand(100:1000))
        h = t3.(rand(100:1000, rand(2:10)))
        r = r0 .+ accumulate(+,sort(h)) |> Vector{t2}

        # outputs
        if t1 <: Integer
            @test_throws MethodError GlobalStructured3DGrid_EqSim(t1, r0, h, rand(2:10,2)...)
            @test_throws MethodError GlobalStructured3DGrid_EqSim(t1, r, rand(2:10,2)...)
        else
            @test GlobalStructured3DGrid_EqSim(t1, r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{t1}
            @test GlobalStructured3DGrid_EqSim(t1, r, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{t1}
        end

        if t2 <: Integer
            if t2 == BigInt
                @test GlobalStructured3DGrid_EqSim(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{BigFloat}
                @test GlobalStructured3DGrid_EqSim(r, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{BigFloat}
            else
                @test GlobalStructured3DGrid_EqSim(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{Float64}
                @test GlobalStructured3DGrid_EqSim(r, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{Float64}
            end
        else
            @test GlobalStructured3DGrid_EqSim(r0, h, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{t2}
            @test GlobalStructured3DGrid_EqSim(r, rand(2:10,2)...) isa GlobalStructured3DGrid_EqSim{t2}
        end
    end

    #::. area and volume tests
    r, h = Float64.(rand(2) .* 1000)
    grid = GlobalStructured3DGrid_EqSim(r, [h], rand(10:20,2)...)
    Nsurface = length(surfacecoords(grid))
    
    @test isapprox(sum(areas(grid)[1:Nsurface]), 2 * pi * r^2; rtol=1e-4)
    @test isapprox(sum(volumes(grid)), 2/3 * pi * ((r+h)^3 - r^3); rtol=1e-4)

    #::. coord2idx
    r, h = rand()*100, sort(rand(10)*100)
    grid = GlobalStructured3DGrid_EqSim(r, h, rand(10:20,2)...)
    Xs = coords(grid)
    for _ in 1:100
        idx = rand(1:length(Xs))
        @test coord2idx(grid, Xs[idx]) == idx
    end

    #::. surfacecoords
    Xssurf = surfacecoords(grid)
    for xs in Xssurf; @test xs.r == r; end

    ########################################################################################
    # GlobalStructured3DGrid_Reduced                                                       #
    ########################################################################################

    #::. type stability and constructor tests
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:100

        # input types
        t1, t2, t3 = rand(types, 3)
        r0 = t2(rand(100:1000))
        h = t3.(rand(100:1000, rand(2:10)))
        r = r0 .+ accumulate(+,sort(h)) |> Vector{t2}

        # outputs
        if t1 <: Integer
            @test_throws MethodError GlobalStructured3DGrid_Reduced(t1, r0, h, rand(2:10))
            @test_throws MethodError GlobalStructured3DGrid_Reduced(t1, r, rand(2:10))
        else
            @test GlobalStructured3DGrid_Reduced(t1, r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced{t1}
            @test GlobalStructured3DGrid_Reduced(t1, r, rand(2:10)) isa GlobalStructured3DGrid_Reduced{t1}
        end

        if t2 <: Integer
            if t2 == BigInt
                @test GlobalStructured3DGrid_Reduced(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced{BigFloat}
                @test GlobalStructured3DGrid_Reduced(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced{BigFloat}
            else
                @test GlobalStructured3DGrid_Reduced(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced{Float64}
                @test GlobalStructured3DGrid_Reduced(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced{Float64}
            end
        else
            @test GlobalStructured3DGrid_Reduced(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced{t2}
            @test GlobalStructured3DGrid_Reduced(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced{t2}
        end
    end

    #::. area and volume tests
    r, h = Float64.(rand(2) .* 1000)
    grid = GlobalStructured3DGrid_Reduced(r, [h], rand(10:20))
    Nsurface = length(surfacecoords(grid))
    
    @test isapprox(sum(areas(grid)[1:Nsurface]), 4 * pi * r^2; rtol=1e-4)
    @test isapprox(sum(volumes(grid)), 4/3 * pi * ((r+h)^3 - r^3); rtol=1e-4)

    #::. coord2idx
    r, h = rand()*100, sort(rand(10)*100)
    grid = GlobalStructured3DGrid_Reduced(r, h, rand(10:20))
    Xs = coords(grid)
    for _ in 1:100
        idx = rand(1:length(Xs))
        @test coord2idx(grid, Xs[idx]) == idx
    end

    #::. surfacecoords
    Xssurf = surfacecoords(grid)
    for xs in Xssurf; @test xs.r == r; end


    ########################################################################################
    # GlobalStructured3DGrid_Reduced_EqSim                                                 #
    ########################################################################################

    #::. type stability and constructor tests
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:100

        # input types
        t1, t2, t3 = rand(types, 3)
        r0 = t2(rand(100:1000))
        h = t3.(rand(100:1000, rand(2:10)))
        r = r0 .+ accumulate(+,sort(h)) |> Vector{t2}

        # outputs
        if t1 <: Integer
            @test_throws MethodError GlobalStructured3DGrid_Reduced_EqSim(t1, r0, h, rand(2:10))
            @test_throws MethodError GlobalStructured3DGrid_Reduced_EqSim(t1, r, rand(2:10))
        else
            @test GlobalStructured3DGrid_Reduced_EqSim(t1, r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{t1}
            @test GlobalStructured3DGrid_Reduced_EqSim(t1, r, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{t1}
        end

        if t2 <: Integer
            if t2 == BigInt
                @test GlobalStructured3DGrid_Reduced_EqSim(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{BigFloat}
                @test GlobalStructured3DGrid_Reduced_EqSim(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{BigFloat}
            else
                @test GlobalStructured3DGrid_Reduced_EqSim(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{Float64}
                @test GlobalStructured3DGrid_Reduced_EqSim(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{Float64}
            end
        else
            @test GlobalStructured3DGrid_Reduced_EqSim(r0, h, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{t2}
            @test GlobalStructured3DGrid_Reduced_EqSim(r, rand(2:10)) isa GlobalStructured3DGrid_Reduced_EqSim{t2}
        end
    end

    #::. area and volume tests
    r, h = Float64.(rand(2) .* 1000)
    grid = GlobalStructured3DGrid_Reduced_EqSim(r, [h], rand(10:20))
    Nsurface = length(surfacecoords(grid))
    
    @test isapprox(sum(areas(grid)[1:Nsurface]), 2 * pi * r^2; rtol=1e-4)
    @test isapprox(sum(volumes(grid)), 4/3 * pi * ((r+h)^3 - r^3) / 2; rtol=1e-4)

    #::. coord2idx
    r, h = rand()*100, sort(rand(10)*100)
    grid = GlobalStructured3DGrid_Reduced_EqSim(r, h, rand(10:20))
    Xs = coords(grid)
    for _ in 1:100
        idx = rand(1:length(Xs))
        @test coord2idx(grid, Xs[idx]) == idx
    end

    #::. surfacecoords
    Xssurf = surfacecoords(grid)
    for xs in Xssurf; @test xs.r == r; end

end

nothing