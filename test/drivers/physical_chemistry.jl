print("TESTING: derivers > physical_chemistry.jl")

@testset verbose=true "physical_chemistry.jl ........." begin

    #::. type-stability
    types = [Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:1000

        # draw random types
        t1, t2, t3, t4 = rand(types, 4)
        t5 = rand(types[4:7])

        # setup inputs
        in1, in2, in3, in4 = rand(4)
        in1 = t1 <: Integer ? round(t1, in1) : t1(in1)
        in2 = t2 <: Integer ? round(t2, in2) : t2(in2)
        in3 = t3 <: Integer ? round(t3, in3) : t3(in3)
        in4 = t4 <: Integer ? round(t4, in4) : t4(in4)

        # prepare output type(s)
        t_out_2 = promote_type(t1, t2)
        if t_out_2 <: Integer; t_out_2 = promote_type(t_out_2, Float64); end
        t_out_3 = promote_type(t1, t2, t3)
        if t_out_3 <: Integer; t_out_3 = promote_type(t_out_3, Float64); end
        t_out_4 = promote_type(t1, t2, t3, t4)
        if t_out_4 <: Integer; t_out_4 = promote_type(t_out_4, Float64); end

        # arrhenius
        @test arrhenius(in1, in2) isa t_out_2
        @test arrhenius(in1, in2, in3) isa t_out_3
        @test arrhenius(t5, in1, in2) isa t5
        @test arrhenius(t5, in1, in2, in3) isa t5

        # diffusion_time
        @test diffusion_time(in1, in2, in3, in4) isa t_out_4
        @test diffusion_time(t5, in1, in2, in3, in4) isa t5
    end


    #::. behaviour


end

println("\rTESTING: drivers > physical_chemistry.jl - DONE")
nothing
