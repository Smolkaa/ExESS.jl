using Test
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "src", "ExESS.jl"))
    using .ExESS
end

#::. begin ExESS package test
@testset verbose=true "ExESS =============================" begin
        
    # base
    @testset verbose=true "BASE ----------------------------" begin

        include(joinpath(@__DIR__, "base", "constants.jl"))
        include(joinpath(@__DIR__, "base", "coordinates.jl"))
        include(joinpath(@__DIR__, "base", "utils.jl"))
    end

    # drivers
    @testset verbose=true "DRIVERS -------------------------" begin
    end

    # grids
    @testset verbose=true "GRIDS ---------------------------" begin
    end

    # misc
    @testset verbose=true "MISC ----------------------------" begin
    end

    # moon
    @testset verbose=true "MOON ----------------------------" begin
    end

    # trajectories
    @testset verbose=true "TRAJECTORIES --------------------" begin
    end

end

nothing