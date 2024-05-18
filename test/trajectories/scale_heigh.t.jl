print("TESTING: trajectories > scale_height.jl")

@testset verbose=true "scale_height.jl ............" begin

    #::. type_stability
    types = [Int32, Int64, BigInt, Float32, Float64, BigFloat]
    for t1 in types, t2 in types

        T = t1 <: Integer ? t1(rand(1:400)) : rand(t1) * 400
        m = t2 <: Integer ? t2(rand(1:10)) : amu2kg(rand(t2)*100)

        t_out = promote_type(t1, t2)
        if t_out <: Integer; t_out = promote_type(t_out, Float64); end

        @test scale_height(T, m) isa t_out
    end

    #::. physical behaviour
    @test all(_ -> (scale_height(0, amu2kg(rand()*10)) == 0), 1:1000)
    @test all(_ -> isinf(scale_height(rand()*400, 0)), 1:1000)
end

println("\rTESTING: trajectories > scale_height.jl - DONE")
nothing
