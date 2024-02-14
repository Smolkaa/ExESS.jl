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
        i, j, k, l = rand(1:8, 4)
        t1, t2, t3, t4, t5 = types[i], types[j], types[k], types[l], types[rand(5:8)]

        # prepare intputs
        in1 = t1 <: Integer ? rand(1:100) |> t1 : rand(t1)
        in2 = t2 <: Integer ? rand(1:100) |> t2 : rand(t2)
        in3 = t3 <: Integer ? rand(1:100) |> t3 : rand(t3)
        in4 = t4 <: Integer ? rand(1:100) |> t4 : rand(t4)

        # prepare velocities vectors, tuples, and custom types
        tp3, tp4 = (in3, in3, in3), (in4, in4, in4)
        v3, v4 = [in3, in3, in3], [in4, in4, in4]
        vl3 = LocalCartesianVelocity(in3, in3, in3)
        vl4 = LocalCartesianVelocity(in4, in4, in4)

        # prepare output types
        t_out_1 = t1 <: Integer ? typeof(promote(one(t1), 1.0)[1]) : t1
        t_out_2 = typeof(promote(in1, in2)[1])
        if t_out_2 <: Integer; t_out_2 = typeof(promote(one(t_out_2), 1.0)[1]); end


        # test all maxwellian distributions
        for d in [
            MBAzimuthDistribution(in1, in2), 
            MBElevationDistribution(in1, in2),
            MBSpeedDistribution(in1, in2),
            MBVelocityDistribution(in1, in2),
            MBFluxAzimuthDistribution(in1, in2),
            MBFluxElevationDistribution(in1, in2),
            MBFluxSpeedDistribution(in1, in2),
            MBFluxVelocityDistribution(in1, in2),]

            # test for type stability of constructor
            @test typeof(d.T) == t_out_2
            @test typeof(d.m) == t_out_2

            # split tests for different input types
            if d isa Union{MBVelocityDistribution, MBFluxVelocityDistribution}
                @test typeof(Base.rand(d))     == Tuple{t_out_2, t_out_2, t_out_2}

                @test typeof(cdf(d, tp3, v4))  == t_out_2
                @test typeof(cdf(d, v3, vl4))  == t_out_2
                @test typeof(cdf(d, vl3, tp4)) == t_out_2
                @test typeof(cdf(t5, d, tp3, v4))  == t5
                @test typeof(cdf(t5, d, v3, vl4))  == t5
                @test typeof(cdf(t5, d, vl3, tp4)) == t5
                
                @test typeof(pdf(d, tp3)) == t_out_2
                @test typeof(pdf(d, v3))  == t_out_2
                @test typeof(pdf(d, vl3)) == t_out_2
                @test typeof(pdf(t5, d, tp3)) == t5
                @test typeof(pdf(t5, d, v3))  == t5
                @test typeof(pdf(t5, d, vl3)) == t5

            else
                @test typeof(Base.rand(d)) == t_out_2

                @test typeof(cdf(d, in3, in4)) == t_out_2
                @test typeof(cdf(t5, d, in3, in4)) == t5

                @test typeof(pdf(d, in3)) == t_out_2
                @test typeof(pdf(t5, d, in3)) == t5

            end
        end
    end



end