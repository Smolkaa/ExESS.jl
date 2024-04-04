#=
    Smallest integer types `Int16` and `Int32` are not included here as it is generally too 
    small for some planetary applications.
=#
@testset "type stability" begin

    # Monte Carlo Test (N=1000)
    types = [Int32, Int64, BigInt, Float16, Float32, Float64, BigFloat]
    for _ in 1:250
        # draw random types
        t1, t2, t3, t4, t5 = rand(types, 5)

        # prepare inputs
        R1 = LUNAR_RADIUS * (1 + 5*rand()) * 1e-6 
        if t1 <: Integer; R1 = ceil(t1, R1); else; R1 = t1(R1); end

        R2 = R1 + LUNAR_RADIUS * rand() * 1e-6
        if t2 <: Integer; R2 = ceil(t2, R2); else; R2 = t2(R2); end

        M  = LUNAR_MASS * 1e-18
        if t3 <: Integer; M = ceil(t3, M); else; M = t3(M); end

        m = amu2kg(100*rand())
        if t4 <: Integer; m = ceil(t4, m); else; m = t4(m); end

        T  = 100*rand()
        if t5 <: Integer; T = ceil(t5, T); else; T = t5(T); end

        # precalculate output type
        tout = promote_type(t1, t2, t3, t4, t5)
        if tout <: Integer; tout = promote_type(tout, Float64); end
        tout2 = rand(types[4:7])

        # test type stability
        @test projection_CHAMBERLAIN1963(R1, R2, M, m, T) isa tout
        @test projection_CHAMBERLAIN1963(R1+R2, R1+R2, M, m, T) isa tout
        @test projection_CHAMBERLAIN1963(tout2, R1, R2, M, m, T) isa tout2
        @test projection_CHAMBERLAIN1963(tout2, R1+R2, R1+R2, M, m, T) isa tout2
    end

end

nothing