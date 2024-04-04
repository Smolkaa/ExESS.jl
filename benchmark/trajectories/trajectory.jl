using BenchmarkTools, Dates, Printf
if !isdefined(Main, :ExESS)
    include(joinpath(@__DIR__, "..", "..", "src", "ExESS.jl"))
    using .ExESS
end

#::. define benchmark
function benchmark_trajectory()

    # open output file
    fid = open(joinpath(@__DIR__, "__trajectory__.md"), "w+")

    # write file name
    write(fid, "# Benchmark: `trajectory` \n\n")

    # write time of last benchmark
    write(fid, "Last updated: $(Dates.now())\n\n")

    # write julia version information
    write(fid, "**Julia Version Information**\n\n")
    write(fid, "```julia\n")
    versioninfo(fid)
    write(fid, "```\n\n")

    # inputs
    x0_gsp = rand(GlobalSphericalPosition, SolarSurfaceDistribution(LUNAR_RADIUS))
    x0_gcp = GlobalCartesianPosition(x0_gsp)
    x0_t = ExESS._get(x0_gcp)
    x0_v = vec(x0_gcp)

    v0_lcv = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
    v0_gcv = GlobalCartesianVelocity(x0_gcp, v0_lcv)
    v0_t = ExESS._get(v0_gcv)
    v0_v = vec(v0_gcv)

    # write inputs to file
    write(fid, "## Preprocessing\n\n")
    write(fid, "```julia\n")
    write(fid, "x0_gsp = rand(GlobalSphericalPosition, SolarSurfaceDistribution(LUNAR_RADIUS))\n")
    write(fid, "x0_gcp = GlobalCartesianPosition(x0_gsp)\n")
    write(fid, "x0_t = ExESS._get(x0_gcp)\n")
    write(fid, "x0_v = vec(x0_gcp)\n")
    write(fid, "\n")
    write(fid, "v0_lcv = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))\n")
    write(fid, "v0_gcv = GlobalCartesianVelocity(x0_gcp, v0_lcv)\n")
    write(fid, "v0_t = ExESS._get(v0_gcv)\n")
    write(fid, "v0_v = vec(v0_gcv)\n")
    write(fid, "```\n\n")

    # benchmarks
    write(fid, "## Benchmarks\n\n")


    write(fid, "### Function: `trajectory` - Tuples\n\n")
    bm = @benchmark trajectory($x0_t, $v0_t, $ddx_gravity) evals=100 seconds=10
    write(fid, "```julia\n")
    write(fid, "trajectory(x0_t, v0_t, ddx_gravity)\n")
    write(fid, "```\n\n")
    @printf(fid, "- Minium time: %.1f ns\n", minimum(bm.times))
    @printf(fid, "- Median time: %.1f ns\n", median(bm.times))
    @printf(fid, "- Memory estimate: %i bytes\n", bm.memory)
    @printf(fid, "- Allocs estimate: %i\n\n", bm.allocs)
    

    write(fid, "### Function: `trajectory` - Tuples w/ `reltol=1e-8`\n\n")
    bm = @benchmark trajectory($x0_t, $v0_t, $ddx_gravity; reltol=1e-8) evals=100 seconds=10
    write(fid, "```julia\n")
    write(fid, "trajectory(x0_t, v0_t, ddx_gravity; reltol=1e-8)\n")
    write(fid, "```\n\n")
    @printf(fid, "- Minium time: %.1f ns\n", minimum(bm.times))
    @printf(fid, "- Median time: %.1f ns\n", median(bm.times))
    @printf(fid, "- Memory estimate: %i bytes\n", bm.memory)
    @printf(fid, "- Allocs estimate: %i\n\n", bm.allocs)


    write(fid, "### Function: `trajectory` - Vectors\n\n")
    bm = @benchmark trajectory($x0_v, $v0_v, $ddx_gravity) evals=100 seconds=10
    write(fid, "```julia\n")
    write(fid, "trajectory(x0_v, v0_v, ddx_gravity)\n")
    write(fid, "```\n\n")
    @printf(fid, "- Minium time: %.1f ns\n", minimum(bm.times))
    @printf(fid, "- Median time: %.1f ns\n", median(bm.times))
    @printf(fid, "- Memory estimate: %i bytes\n", bm.memory)
    @printf(fid, "- Allocs estimate: %i\n\n", bm.allocs)


    write(fid, "### Function: `trajectory` - Custom Types\n\n")
    bm = @benchmark trajectory($x0_gcp, $v0_gcv, $ddx_gravity) evals=100 seconds=10
    write(fid, "```julia\n")
    write(fid, "trajectory(x0_gcp, v0_gcv, ddx_gravity)\n")
    write(fid, "```\n\n")
    @printf(fid, "- Minium time: %.1f ns\n", minimum(bm.times))
    @printf(fid, "- Median time: %.1f ns\n", median(bm.times))
    @printf(fid, "- Memory estimate: %i bytes\n", bm.memory)
    @printf(fid, "- Allocs estimate: %i\n\n", bm.allocs)


    write(fid, "### Function: `trajectory` - Custom Types w/ conversion\n\n")
    bm = @benchmark trajectory($x0_gsp, $v0_lcv, $ddx_gravity) evals=100 seconds=10
    write(fid, "```julia\n")
    write(fid, "trajectory(x0_gsp, v0_lcv, ddx_gravity)\n")
    write(fid, "```\n\n")
    @printf(fid, "- Minium time: %.1f ns\n", minimum(bm.times))
    @printf(fid, "- Median time: %.1f ns\n", median(bm.times))
    @printf(fid, "- Memory estimate: %i bytes\n", bm.memory)
    @printf(fid, "- Allocs estimate: %i\n\n", bm.allocs)


    close(fid)
    nothing
end
benchmark_trajectory()