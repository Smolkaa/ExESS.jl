using BenchmarkTools, Dates, Printf, LinearAlgebra
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "..", "src", "ExESS.jl"))
    using .ExESS
end
using DifferentialEquations


function benchmark_ode_algorithms(alg)

    #
    mintime = Float32[]
    for logv in range(0, 3, length=7)
        v = 10^logv

        x0 = rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS))
        v0 = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
        v0 = normalize(v0) * v
        bm = @benchmark trajectory($x0, $v0; alg=$alg) evals=100 seconds=20
        push!(mintime, minimum(bm.times))
    end

    #
    return mintime
end
