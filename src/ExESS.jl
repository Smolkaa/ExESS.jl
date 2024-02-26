module ExESS


#::. loading packages
using DelimitedFiles: readdlm
using DifferentialEquations: ODESolution, SecondOrderODEProblem, solve, terminate!, Tsit5, 
                             VectorContinuousCallback
using Interpolations: BSpline, Flat, interpolate, Linear, linear_interpolation
using LinearAlgebra
using NearestNeighbors: KDTree, knn
using SpecialFunctions: gamma, gamma_inc, erf
using StaticArrays: SA
using Statistics


#::. BASE
include(joinpath(@__DIR__, "base", "constants.jl"))
include(joinpath(@__DIR__, "base", "coordinates.jl"))
include(joinpath(@__DIR__, "base", "statistics.jl"))
include(joinpath(@__DIR__, "base", "utils.jl"))

#::. DRIVERS
include(joinpath(@__DIR__, "drivers", "maxwellian.jl"))
include(joinpath(@__DIR__, "drivers", "surface_distributions.jl"))
include(joinpath(@__DIR__, "drivers", "thermal_sorption.jl"))


#::. GRIDS
include(joinpath(@__DIR__, "grids", "grids.jl"))
include(joinpath(@__DIR__, "grids", "spherical2d.jl"))
include(joinpath(@__DIR__, "grids", "spherical2d_spiral.jl"))
include(joinpath(@__DIR__, "grids", "spherical2d_HEALPix.jl"))
include(joinpath(@__DIR__, "grids", "spherical3d.jl"))
include(joinpath(@__DIR__, "grids", "cartesian3d.jl"))


#::. TRAJECTORY
include(joinpath(@__DIR__, "trajectories", "trajectory.jl"))
include(joinpath(@__DIR__, "trajectories", "escape_velocity.jl"))
include(joinpath(@__DIR__, "trajectories", "scale_height.jl"))

#::. CERES
include(joinpath(@__DIR__, "ceres", "constants.jl"))

#::. MERCURY
include(joinpath(@__DIR__, "mercury", "constants.jl"))

#::. MOON
include(joinpath(@__DIR__, "moon", "constants.jl"))
include(joinpath(@__DIR__, "moon", "lunar_surface_temperatures_DIVINER.jl"))
include(joinpath(@__DIR__, "moon", "lunar_surface_temperatures_BUTLER1997.jl"))
include(joinpath(@__DIR__, "moon", "lunar_surface_temperatures_HURLEY2015.jl"))

#::. MISC
include(joinpath(@__DIR__, "misc", "projection_CHAMBERLAIN1963.jl"))
include(joinpath(@__DIR__, "misc", "solar_angle.jl"))




"""
    ExESS
"""
ExESS

end # module
