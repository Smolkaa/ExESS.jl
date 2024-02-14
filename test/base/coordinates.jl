@testset verbose=true "coordinates.jl ................" begin

    include(joinpath(@__DIR__, "coordinates", "type_stability_constructors.jl"))
    include(joinpath(@__DIR__, "coordinates", "type_stability_subtypes.jl"))

    # TODO: general extensions: norm, dot, cross, eltype, etc.
    include(joinpath(@__DIR__, "coordinates", "base_extensions_position.jl"))
    include(joinpath(@__DIR__, "coordinates", "base_extensions_velocity.jl"))

    # TODO: conversion tests >> back and forth >> error estimation
    # TODO: utility tests (incl. type stability) >> azimuth, elevation, speed, _get, etc.
end

nothing