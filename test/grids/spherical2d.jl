@testset verbose=true "spherical2d.jl ................" begin

    include(joinpath(@__DIR__, "spherical2d", "type_stability.jl")) 
    include(joinpath(@__DIR__, "spherical2d", "behaviour.jl")) 

end

nothing