@testset verbose=true "trajectory.jl .............." begin
    print("TESTING: trajectories > trajectory.jl")

    #
    include(joinpath(@__DIR__, "trajectory", "type_stability.jl"))
    include(joinpath(@__DIR__, "trajectory", "behaviour.jl"))


    println("\rTESTING: trajectories > trajectory.jl - DONE")
end

nothing
