#=
    The type stability of the maxwellian distributions is usually goverend by the subtype
    of the individual distribution type. In the `base/statistics.jl` file, an additional
    method is defined that takes an additional input type as first argument, which
    overwrites the type of the output. This is tested here as well.
=#

@testset "type stability" begin

    # Monte Carlo tests (N=100)
    types = [Int16, Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:100

        # draw random types
        t1, t2, t3 = rand(types, 3)
        t4 = rand(types[5:8])
        t_out_1 = t1 <: Integer ? promote_type(t1, Float64) : t1

        # prepare inputs
        r = t1 <: Integer ? t1(rand(1:100)) : t1(rand() * 100)
        l_t = t2 <: Integer ? t2.((1, rand(-1:1), rand(-1:1))) : 
                              t2.((1, rand() * pi - pi/2, rand() * pi - pi/2))
        u_t = t3 <: Integer ? t3.((1, rand(-1:1), rand(-1:1))) : 
                              t3.((1, rand() * pi - pi/2, rand() * pi - pi/2))        
        l_v, u_v = [l_t[1], l_t[2], l_t[3]], [u_t[1], u_t[2], u_t[3]]
        l_s, u_s = GlobalSphericalPosition(l_t), GlobalSphericalPosition(u_t)

        # test all surface distributions
        for d in [EqualSurfaceDistribution(r), SolarSurfaceDistribution(r)]

            # Base.rand
            @test Base.rand(d) isa Tuple{t_out_1, t_out_1, t_out_1}

            # cdf
            @test cdf(d, l_t, u_v)  isa t_out_1
            @test cdf(d, l_v, u_s)  isa t_out_1
            @test cdf(d, l_s, u_t)  isa t_out_1
            @test cdf(t4, d, l_t, u_v)  isa t4
            @test cdf(t4, d, l_v, u_s)  isa t4
            @test cdf(t4, d, l_s, u_t)  isa t4

            # pdf
            @test pdf(d, l_t)   isa t_out_1
            @test pdf(d, l_v)   isa t_out_1
            @test pdf(d, l_s)   isa t_out_1
            @test pdf(t4, d, u_t)   isa t4
            @test pdf(t4, d, u_v)   isa t4
            @test pdf(t4, d, u_s)   isa t4

        end
    end



end