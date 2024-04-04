@testset "behaviour - partitions" begin

    # Page 911 - partition function of bound orbits
    # Testing the sum of the ballistic and satellite partition functions
    # 
    # Note that the base radius `R0` is divided by ten to avoid R < R0.
    R0, m, M, T = LUNAR_RADIUS/10, amu2kg(100*rand()), LUNAR_MASS, 400*rand()
    LMB = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * R0)

    lmb = vcat(0.1:0.1:1.0, 1.5:0.5:5.0, 6.0:1.0:15.0)
    P = [
        0.0224107,  # 0.1
        0.0597575,  # 0.2
        0.1035676,  # 0.3
        0.1505329,  # 0.4
        0.1987480,  # 0.5
        0.2469956,  # 0.6
        0.2944652,  # 0.7
        0.3406101,  # 0.8
        0.3850650,  # 0.9
        0.4275932,  # 1.0
        
        0.6083748,  # 1.5
        0.7385358,  # 2.0
        0.8282028,  # 2.5
        0.8883897,  # 3.0
        0.9281022,  # 3.5
        0.9539882,  # 4.0
        0.9707090,  # 4.5
        0.9814339,  # 5.0

        0.9926172,  # 6.0
        0.9970952,  # 7.0
        0.9988663,  # 8.0
        0.9995605,  # 9.0
        0.9998306,  # 10.0
        0.9999351,  # 11.0
        0.9999754,  # 12.0
        0.9999908,  # 13.0
        0.9999967,  # 14.0
        0.9999990,  # 15.0
    ]

    for (i, l) in enumerate(lmb)
        R = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * l)
        @test isapprox(projection_CHAMBERLAIN1963(R0, R, M, m, T; esc=false) * exp(LMB - l),  
                       P[i]; atol=1e-6)
    end


    # Page 912 - partition function of ballistic and escaping orbits
    # Testing ballistic and escaping partition fraction individually
    LMB1 = [1.0, 2.0, 3.0, 4.0, 5.0, 7.5]
    LMB2 = vcat(0.1:0.1:0.5)
    PBAL = [
        # 1.0       # 2.0       # 3.0       # 4.0       # 5.0       # 7.5
        .0031541    .0016478    .0011154    .0008430    .0006775    .0004545;   # 0.1
        .0158253    .0086060    .0059134    .0045046    .0036380    .0024567;   # 0.2
        .0388161    .0218891    .0152638    .0117200    .0095123    .0064675;   # 0.3
        .0711837    .0414704    .0293366    .0227057    .0185221    .0126828;   # 0.4
        .1115230    .0668653    .0479653    .0374208    .0306831    .0211637;   # 0.5
    ]
    PESC = [
        # 1.0       # 2.0       # 3.0       # 4.0       # 5.0       # 7.5
        .0054313    .0021737    .0013295    .0009523    .0007404    .0004745;   # 0.1
        .0182501    .0071668    .0043529    .0031064    .0024094    .0015391;   # 0.2
        .0355262    .0136826    .0082590    .0058754    .0045483    .0028978;   # 0.3
        .0556409    .0209618    .0125751    .0089191    .0068921    .0043806;   # 0.4
        .0777011    .0285147    .0169971    .0120195    .0092718    .0058796;   # 0.5
    ]

    for (i, l1) in enumerate(LMB1)
        for (j, l2) in enumerate(LMB2)
            R1 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * l1)
            R2 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * l2)
            @test isapprox(
                projection_CHAMBERLAIN1963(R1, R2, M, m, T; esc=false, sat=false)*exp(l1 - l2), 
                PBAL[j, i]; atol=1e-6)
            @test isapprox(
                projection_CHAMBERLAIN1963(R1, R2, M, m, T; bal=false, sat=false)*exp(l1 - l2), 
                PESC[j, i]; atol=1e-6)
        end
    end

end

nothing