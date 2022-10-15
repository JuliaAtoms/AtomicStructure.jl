@testset "Potentials" begin
    @testset "Point charges" begin
        @test charge(pc"H") == 1
        @test charge(pc"He") == 2
        @test charge(pc"Xe") == 54
        @test effective_charge(pc"Xe") == 54

        @test string(pc"He") == "Z = 2 [helium]"
        @test string(PointCharge(3.4)) == "Z = 3.4"

        @test pc"H"(o"1s", [1,2]) == -[1, 0.5]

        @test ground_state(pc"H") == c"1s"
        @test ground_state(pc"He") == c"1s2"

        @test islocal(pc"H")
        @test isspinlocal(pc"H")

        @test pc"H"+3 == pc"Be"
        @test pc"Be"-2 == pc"He"
    end

    @testset "Yukawa" begin
        V = Yukawa(1.5, 1.0, 2.0)
        @test charge(V) == 2.25
        @test effective_charge(V) == 2.25
        @test ground_state(V) == c"1s"

        @test islocal(V)
        @test isspinlocal(V)

        T = typeof(V.g)
        @test string(V) == "Yukawa{$T}(1.500, 1.000, 2.000): -2.250*exp(-2.000*r)/r"

        @test V(o"1s", [1,2]) ≈ [-0.3045043872823786, -0.02060509374982595]
    end
end
