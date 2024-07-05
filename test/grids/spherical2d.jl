@testset verbose=true "spherical2d.jl ................" begin
    print("TESTING: grids > spherical2d.jl")

    #::. tests
    include(joinpath(@__DIR__, "spherical2d", "type_stability.jl"))
    include(joinpath(@__DIR__, "spherical2d", "behaviour.jl"))


    println("\rTESTING: grids > spherical2d.jl - DONE")
end

nothing
