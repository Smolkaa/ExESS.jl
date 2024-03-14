@testset verbose=true "trajectory.jl .............." begin

    include(joinpath(@__DIR__, "trajectory", "type_stability.jl"))
    include(joinpath(@__DIR__, "trajectory", "behaviour.jl"))

end

nothing