@testset "behaviour - typical speeds" begin

    #::. typical speeds of velocity distributions (physical behavior tests)
    k_B, T, m = BOLTZMANN_CONSTANT, 400*rand(), amu2kg(100*rand())

    v_mbvd_p, v_mbfvd_p = sqrt(2*k_B*T/m), sqrt(3*k_B*T/m)
    PDF_mbvd = [pdf(MBSpeedDistribution(T, m), v) for v in 1:100_000]
    PDF_mbfvd = [pdf(MBFluxSpeedDistribution(T, m), v) for v in 1:100_000]
    @test isapprox(v_mbvd_p, findmax(PDF_mbvd)[2]; rtol=1e-2)
    # @test isapprox(v_mbvd_p, speed(MaxwellBoltzmannVelocityDistribution(T, m), :prob); rtol=1e-2)
    @test isapprox(v_mbfvd_p, findmax(PDF_mbfvd)[2]; rtol=1e-2)
    # @test isapprox(v_mbfvd_p, speed(MaxwellBoltzmannFluxVelocityDistribution(T, m), :prob); rtol=1e-2)

    v_mbvd_mean, v_mbfvd_mean = sqrt(8*k_B*T/(pi*m)), sqrt(9*pi*k_B*T/(8*m))
    v_mbvd_rms, v_mbfvd_rms = sqrt(3*k_B*T/m), sqrt(4*k_B*T/m)
    V_mbvd = [rand(MBSpeedDistribution(T, m)) for _ in 1:100_000]
    V_mbfvd = [rand(MBFluxSpeedDistribution(T, m)) for _ in 1:100_000]
    @test isapprox(v_mbvd_mean, mean(V_mbvd); rtol=1e-2)
    # @test isapprox(v_mbvd_mean, speed(MaxwellBoltzmannVelocityDistribution(T, m), :mean); rtol=1e-2)
    @test isapprox(v_mbfvd_mean, mean(V_mbfvd); rtol=1e-2)
    # @test isapprox(v_mbfvd_mean, speed(MaxwellBoltzmannFluxVelocityDistribution(T, m), :mean); rtol=1e-2)
    @test isapprox(v_mbvd_rms, sqrt(mean(V_mbvd.^2)); rtol=1e-2)
    # @test isapprox(v_mbvd_rms, speed(MaxwellBoltzmannVelocityDistribution(T, m), :rms); rtol=1e-2)
    @test isapprox(v_mbfvd_rms, sqrt(mean(V_mbfvd.^2)); rtol=1e-2)
    # @test isapprox(v_mbfvd_rms, speed(MaxwellBoltzmannFluxVelocityDistribution(T, m), :rms); rtol=1e-2)

    #::. typical angles
    el_mb_mean, el_mbf_mean = (pi-2)/2, pi/4
    el_mb = abs.([rand(MBElevationDistribution(T, m)) for _ in 1:100_000])
    el_mbf = [rand(MBFluxElevationDistribution(T, m)) for _ in 1:100_000]
    @test isapprox(el_mb_mean, mean(el_mb); rtol=1e-2)
    @test isapprox(el_mbf_mean, mean(el_mbf); rtol=1e-2)   

end