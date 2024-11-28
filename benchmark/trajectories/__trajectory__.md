# Benchmark: `trajectory` 

Last updated: 2024-11-28T09:12:08.356

**Julia Version Information**

```julia
Julia Version 1.11.1
Commit 8f5b7ca12a (2024-10-16 10:53 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: 16 × 13th Gen Intel(R) Core(TM) i7-1360P
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, goldmont)
Threads: 1 default, 0 interactive, 1 GC (on 16 virtual cores)
```

## Preprocessing

```julia
x0_gsp = rand(GlobalSphericalPosition, SolarSurfaceDistribution(LUNAR_RADIUS))
x0_gcp = GlobalCartesianPosition(x0_gsp)
x0_t = Tuple(x0_gcp)
x0_v = vec(x0_gcp)

v0_lcv = rand(LocalCartesianVelocity, MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
v0_gcv = GlobalCartesianVelocity(x0_gcp, v0_lcv)
v0_t = Tuple(v0_gcv)
v0_v = vec(v0_gcv)
```

## Benchmarks

### Function: `trajectory` - Tuples

```julia
trajectory(x0_t, v0_t, ddx_gravity)
```

- Minium time: 88994.0 ns
- Median time: 104287.0 ns
- Memory estimate: 78656 bytes
- Allocs estimate: 1625

### Function: `trajectory` - Tuples w/ `reltol=1e-8`

```julia
trajectory(x0_t, v0_t, ddx_gravity; reltol=1e-8)
```

- Minium time: 107582.0 ns
- Median time: 134363.0 ns
- Memory estimate: 130752 bytes
- Allocs estimate: 2697

### Function: `trajectory` - Vectors

```julia
trajectory(x0_v, v0_v, ddx_gravity)
```

- Minium time: 90028.0 ns
- Median time: 102799.5 ns
- Memory estimate: 78880 bytes
- Allocs estimate: 1635

### Function: `trajectory` - Custom Types

```julia
trajectory(x0_gcp, v0_gcv, ddx_gravity)
```

- Minium time: 89812.0 ns
- Median time: 113712.0 ns
- Memory estimate: 78656 bytes
- Allocs estimate: 1625

### Function: `trajectory` - Custom Types w/ conversion

```julia
trajectory(x0_gsp, v0_lcv, ddx_gravity)
```

- Minium time: 89657.0 ns
- Median time: 111642.0 ns
- Memory estimate: 78272 bytes
- Allocs estimate: 1617

