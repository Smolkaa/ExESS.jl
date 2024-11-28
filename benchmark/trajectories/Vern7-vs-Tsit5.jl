using BenchmarkTools, Dates, Printf
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "..", "src", "ExESS.jl"))
    using .ExESS
end
using DifferentialEquations


function b()

    #
    mintime = Float32[]
    for logv in range(0, 3, length=25)
        v = 10^logv
        bm = @benchmark __benchmark_trajectory__(Tsit5(), $v) evals=100 seconds=10
        push!(mintime, minimum(bm.times))
    end

    #
    return mintime
end


function __benchmark_trajectory__(alg, v)
    x0 = rand(GlobalSphericalPosition, EqualSurfaceDistribution(LUNAR_RADIUS))
    v0 = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
    v0 = normalize(v0) * v
    _ = trajectory(x0, v0; alg=alg)
    return nothing
end
