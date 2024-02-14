using LinearAlgebra, Statistics

@testset verbose=true "maxwellian.jl ................." begin

    include(joinpath(@__DIR__, "maxwellian", "type_stability.jl")) 
    include(joinpath(@__DIR__, "maxwellian", "behaviour.jl")) 
    include(joinpath(@__DIR__, "maxwellian", "behaviour_typical_speeds.jl")) 

end

nothing
