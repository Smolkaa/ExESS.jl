# Benchmark: `trajectory` 

Last updated: 2024-03-09T21:10:20.701

**Julia Version Information**

```julia
Julia Version 1.10.0
Commit 3120989f39 (2023-12-25 18:01 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: 16 × 13th Gen Intel(R) Core(TM) i7-1360P
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, goldmont)
  Threads: 1 on 16 virtual cores
```

## Preprocessing

```julia
x0_gsp = rand(SolarSurfaceDistribution(LUNAR_RADIUS))
x0_gcp = GlobalCartesianPosition(x0_gsp)
x0_t = ExESS._get(x0_gcp)
x0_v = vec(x0_gcp)

v0_lcv = rand(MBFluxVelocityDistribution(100 + rand()*300, amu2kg(rand()*20)))
v0_gcv = GlobalCartesianVelocity(x0_gcp, v0_lcv)
v0_t = ExESS._get(v0_gcv)
v0_v = vec(v0_gcv)
```

## Benchmarks

### Function: `trajectory` - Tuples

```julia
trajectory(x0_t, v0_t, ddx_gravity)
```

- Minium time: 56979.0 ns
- Median time: 62296.0 ns
- Memory estimate: 79648 bytes
- Allocs estimate: 1164

### Function: `trajectory` - Tuples w/ `reltol=1e-8`

```julia
trajectory(x0_t, v0_t, ddx_gravity; reltol=1e-8)
```

- Minium time: 84599.0 ns
- Median time: 111414.0 ns
- Memory estimate: 140272 bytes
- Allocs estimate: 1863

### Function: `trajectory` - Vectors

```julia
trajectory(x0_v, v0_v, ddx_gravity)
```

- Minium time: 58859.0 ns
- Median time: 107375.0 ns
- Memory estimate: 79872 bytes
- Allocs estimate: 1174

### Function: `trajectory` - Custom Types

```julia
trajectory(x0_gcp, v0_gcv, ddx_gravity)
```

- Minium time: 98356.0 ns
- Median time: 106940.0 ns
- Memory estimate: 79648 bytes
- Allocs estimate: 1164

### Function: `trajectory` - Custom Types w/ conversion

```julia
trajectory(x0_gsp, v0_lcv, ddx_gravity)
```

- Minium time: 99263.0 ns
- Median time: 110141.5 ns
- Memory estimate: 79648 bytes
- Allocs estimate: 1164

