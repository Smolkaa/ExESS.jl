@testset verbose=true "temperatures.jl" begin
    print("TESTING: moon > temperatures.jl")


    #::. lunar_surface_temperatures_DIVINER
    include(joinpath(@__DIR__, "temperatures", "lunar_surface_temperatures_DIVINER.jl"))
    include(joinpath(@__DIR__, "temperatures", "lunar_surface_temperatures_HURLEY2015.jl"))
    include(joinpath(@__DIR__, "temperatures", "lunar_surface_temperatures_BUTLER1997.jl"))



    println("\rTESTING: moon > temperatures.jl - DONE")
end

nothing
