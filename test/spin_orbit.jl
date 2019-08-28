@testset "Spin--orbit interaction" begin
    @testset "Quantum numbers" begin
        orbitals = sos"1[s]"
        @test Atoms.mⱼ(orbitals[1]) == half(1)
        @test Atoms.mⱼ(orbitals[2]) == half(-1)
        
        orbitals = sos"2[p]"
        @test Atoms.mⱼ(orbitals[1]) == half(-1)
        @test Atoms.mⱼ(orbitals[2]) == half(-3)
        @test Atoms.mⱼ(orbitals[3]) == half(1)
        @test Atoms.mⱼ(orbitals[4]) == half(-1)
        @test Atoms.mⱼ(orbitals[5]) == half(3)
        @test Atoms.mⱼ(orbitals[6]) == half(1)
    end
end
