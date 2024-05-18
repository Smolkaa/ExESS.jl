using LinearAlgebra, Statistics

print("TESTING: drivers > surface_distributions.jl")

@testset verbose=true "surface_distributions.jl ......" begin

    include(joinpath(@__DIR__, "surface_distributions", "type_stability.jl"))
    include(joinpath(@__DIR__, "surface_distributions", "behaviour.jl"))

end

println("\rTESTING: drivers > surface_distributions.jl - DONE")
nothing
