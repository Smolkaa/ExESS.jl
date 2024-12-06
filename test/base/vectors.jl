@testset verbose=true "vectors.jl ..............................." begin
    print("TESTING: base > vectors.jl")

    # type stabilities
    include(joinpath(@__DIR__, "vectors", "type_stability_constructors.jl"))
    include(joinpath(@__DIR__, "vectors", "type_stability_subtypes.jl"))

    # TODO: general extensions: norm, dot, cross, eltype, etc.
    include(joinpath(@__DIR__, "vectors", "base_extensions_position.jl"))
    include(joinpath(@__DIR__, "vectors", "base_extensions_velocity.jl"))

    # TODO: conversion tests >> back and forth >> error estimation
    # TODO: utility tests (incl. type stability) >> azimuth, elevation, speed, _get, etc.
    include(joinpath(@__DIR__, "vectors", "utility_functions.jl"))


    println("\rTESTING: base > vectors.jl - DONE")
end

nothing
