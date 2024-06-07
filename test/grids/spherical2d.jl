print("TESTING: grids > spherical2d.jl")

@testset verbose=true "spherical2d.jl ................" begin

    include(joinpath(@__DIR__, "spherical2d", "type_stability.jl"))
    include(joinpath(@__DIR__, "spherical2d", "behaviour.jl"))

end

println("\rTESTING: grids > spherical2d.jl - DONE")
nothing
