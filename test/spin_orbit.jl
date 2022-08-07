@testset "Spin--orbit interaction" begin
    @testset "Quantum numbers" begin
        orbitals = sos"1[s]"
        @test AtomicStructure.mⱼ(orbitals[1]) == half(1)
        @test AtomicStructure.mⱼ(orbitals[2]) == half(-1)
        
        orbitals = sos"2[p]"
        @test AtomicStructure.mⱼ(orbitals[1]) == half(-1)
        @test AtomicStructure.mⱼ(orbitals[2]) == half(-3)
        @test AtomicStructure.mⱼ(orbitals[3]) == half(1)
        @test AtomicStructure.mⱼ(orbitals[4]) == half(-1)
        @test AtomicStructure.mⱼ(orbitals[5]) == half(3)
        @test AtomicStructure.mⱼ(orbitals[6]) == half(1)
    end
end
