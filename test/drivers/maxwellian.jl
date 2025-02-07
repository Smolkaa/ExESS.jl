using LinearAlgebra, Statistics

print("TESTING: drivers > maxwellian.jl")

@testset verbose=true "maxwellian.jl ............................" begin

    include(joinpath(@__DIR__, "maxwellian", "type_stability.jl"))
    include(joinpath(@__DIR__, "maxwellian", "behaviour.jl"))
    include(joinpath(@__DIR__, "maxwellian", "behaviour_typical_speeds.jl"))

end

println("\rTESTING: drivers > maxwellian.jl - DONE")
nothing
