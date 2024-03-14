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
        # TODO: statistics test?
        include(joinpath(@__DIR__, "base", "utils.jl"))
    end

    # drivers
    @testset verbose=true "DRIVERS -------------------------" begin
        include(joinpath(@__DIR__, "drivers", "maxwellian.jl"))
        include(joinpath(@__DIR__, "drivers", "surface_distributions.jl"))
    end

    # grids
    @testset verbose=true "GRIDS ---------------------------" begin
        include(joinpath(@__DIR__, "grids", "spherical2d.jl"))
    end

    # misc
    @testset verbose=true "MISC ----------------------------" begin
    end

    # ceres
    @testset verbose=true "CERES ---------------------------" begin
        include(joinpath(@__DIR__, "ceres", "constants.jl"))
    end

    # mercury
    @testset verbose=true "MERCURY -------------------------" begin
        include(joinpath(@__DIR__, "mercury", "constants.jl"))
    end

    # moon
    @testset verbose=true "MOON ----------------------------" begin
        include(joinpath(@__DIR__, "moon", "constants.jl"))   
    end
    # trajectories
    @testset verbose=true "TRAJECTORIES --------------------" begin
        include(joinpath(@__DIR__, "trajectories", "trajectory.jl"))
    end

end

nothing