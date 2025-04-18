using Test
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "src", "ExESS.jl"))
    using .ExESS
end

#::. begin ExESS package test
@testset verbose=true "ExESS ========================================" begin

    # base
    @testset verbose=true "BASE ---------------------------------------" begin

        include(joinpath(@__DIR__, "base", "constants.jl"))
        include(joinpath(@__DIR__, "base", "vectors.jl"))
        # TODO: statistics test?
        include(joinpath(@__DIR__, "base", "utils.jl"))
    end

    # drivers
    @testset verbose=true "DRIVERS ------------------------------------" begin
        include(joinpath(@__DIR__, "drivers", "maxwellian.jl"))
        include(joinpath(@__DIR__, "drivers", "surface_distributions.jl"))
    end

    # grids
    @testset verbose=true "GRIDS --------------------------------------" begin
        include(joinpath(@__DIR__, "grids", "spherical2d.jl"))
    end

    # misc
    @testset verbose=true "MISC ---------------------------------------" begin
        include(joinpath(@__DIR__, "misc", "projection_CHAMBERLAIN1963.jl"))
        include(joinpath(@__DIR__, "misc", "solar_incidence_angle.jl"))
    end

    # ceres
    @testset verbose=true "CERES --------------------------------------" begin
        include(joinpath(@__DIR__, "ceres", "constants.jl"))
    end

    # mercury
    @testset verbose=true "MERCURY ------------------------------------" begin
        include(joinpath(@__DIR__, "mercury", "constants.jl"))
    end

    # moon
    @testset verbose=true "MOON ---------------------------------------" begin
        include(joinpath(@__DIR__, "moon", "constants.jl"))
        include(joinpath(@__DIR__, "moon", "temperatures.jl"))
    end
    # trajectories
    @testset verbose=true "TRAJECTORIES -------------------------------" begin
        include(joinpath(@__DIR__, "trajectories", "landing_position.jl"))
        include(joinpath(@__DIR__, "trajectories", "trajectory.jl"))
        include(joinpath(@__DIR__, "trajectories", "scale_height.jl"))
    end

    println("")
end

nothing
