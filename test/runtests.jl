using Test
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "src", "ExESS.jl"))
    using .ExESS
end

#::. begin ExESS package test
@testset verbose=true "ExESS ===========================" begin


    include(joinpath(@__DIR__, "constants.jl"))
    include(joinpath(@__DIR__, "coordinates.jl"))



    
#::. base
@testset verbose=true "BASE --------------------------" begin
    # include(joinpath(@__DIR__, "tests", "base", "cs.jl"))
    # include(joinpath(@__DIR__, "tests", "base", "distributions.jl"))
    # include(joinpath(@__DIR__, "tests", "base", "utils.jl"))
end

#::. grids
@testset verbose=true "GRIDS -------------------------" begin
    # include(joinpath(@__DIR__, "tests", "grids", "global_structured_2d_grids.jl"))
    # include(joinpath(@__DIR__, "tests", "grids", "global_structured_3d_grids.jl"))
    # include(joinpath(@__DIR__, "tests", "grids", "local_structured_3d_grids.jl"))
end

#::. exospheres
@testset verbose=true "EXOSPHERES --------------------" begin
    # include(joinpath(@__DIR__, "tests", "exospheres", "landing_position.jl"))
    # include(joinpath(@__DIR__, "tests", "exospheres", "trajectory.jl"))
end

#::. surfaces
@testset verbose=true "SURFACES ----------------------" begin
    # include(joinpath(@__DIR__, "tests", "surfaces", "temperatures.jl"))

end

end

nothing