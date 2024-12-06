print("TESTING: mercury > constants.jl")

@testset verbose=true "constants.jl ............................." begin
    @test isdefined(ExESS, :MERCURY_RADIUS)
    @test isdefined(ExESS, :MERCURY_MASS)
end

println("\rTESTING: mercury > constants.jl - DONE")
nothing
