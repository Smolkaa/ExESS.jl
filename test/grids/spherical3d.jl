@testset verbose=true "spherical3d.jl ................" begin
    print("TESTING: grids > spherical3d.jl")

    #::. tests
    # include(joinpath(@__DIR__, "spherical3d", "type_stability.jl"))
    include(joinpath(@__DIR__, "spherical3d", "behaviour.jl"))


    println("\rTESTING: grids > spherical3d.jl - DONE")
end

nothing
