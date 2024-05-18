print("TESTING: trajectories > trajectory.jl")

@testset verbose=true "trajectory.jl .............." begin

    include(joinpath(@__DIR__, "trajectory", "type_stability.jl"))
    include(joinpath(@__DIR__, "trajectory", "behaviour.jl"))

end

print("\rTESTING: trajectories > trajectory.jl - DONE")
nothing
