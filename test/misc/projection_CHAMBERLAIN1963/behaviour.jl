@testset "behaviour - general" begin

    # r1 < r2: projection smaller than one
    @test all(_ -> (
        projection_CHAMBERLAIN1963(
            LUNAR_RADIUS, 
            LUNAR_RADIUS*(1+rand()),
            LUNAR_MASS,
            amu2kg(100*rand()),
            400*rand(),
        ) < 1.0
    ), 1:1000)

    # r1 > r2: projection bigger than one
    @test all(_ -> (
        projection_CHAMBERLAIN1963(
            LUNAR_RADIUS*(1+rand()),
            LUNAR_RADIUS, 
            LUNAR_MASS,
            amu2kg(100*rand()),
            400*rand(),
        ) > 1.0
    ), 1:1000)

    # r1 == r2: projection equals one
    r = rand()
    @test all(_ -> (
        projection_CHAMBERLAIN1963(
            LUNAR_RADIUS*(1+r),
            LUNAR_RADIUS*(1+r), 
            LUNAR_MASS,
            amu2kg(100*rand()),
            400*rand(),
        ) == 1.0
    ), 1:1000)

end

nothing