print("TESTING: misc > projection_CHAMBERLAIN1963.jl")

@testset verbose=true "projection_CHAMBERLAIN1963.jl ............" begin

    include(joinpath(@__DIR__, "projection_CHAMBERLAIN1963", "type_stability.jl"))
    include(joinpath(@__DIR__, "projection_CHAMBERLAIN1963", "behaviour.jl"))
    include(joinpath(@__DIR__, "projection_CHAMBERLAIN1963", "behaviour_partitions.jl"))

end

println("\rTESTING: misc > projection_CHAMBERLAIN1963.jl - DONE")
nothing
