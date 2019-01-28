using Atoms
using AtomicLevels
using AngularMomentumAlgebra
using SCF
# Needed as /an/ implementation of an AbstractQuasiMatrix, not the
# only possibility:
using FiniteDifferencesQuasi
using Test

@testset "Structure setup" begin
    rₘₐₓ = 300
    ρ = 0.25
    N = ceil(Int, rₘₐₓ/ρ + 1/2)

    @testset "Hydrogenic orbitals" begin
        @testset "Hydrogen" begin
            nucleus = pc"H"
            R = RadialDifferences(N, ρ)
            for method in [:arnoldi_shift_invert, :arnoldi, :eigen]
                atom = Atom(R, csfs([c"1s", c"2s", c"2p", c"3s", c"3p", c"3d"]),
                            nucleus, method=method, verbosity=4)
            end
        end
        @testset "Helium" begin
            nucleus = pc"He"
            R=RadialDifferences(N, ρ, charge(nucleus))
            @testset "CSFs" begin
                atom = Atom(R, csfs([c"1s2", c"1s 2s", c"1s 2p", c"3s 3p", c"4s 3d", c"5s 5d"]),
                            nucleus, verbosity=4)
            end
            @testset "Spin-configurations" begin
                atom = Atom(R, spin_configurations([c"1s2", c"1s 2s", c"1s 2p", c"3s 3p", c"4s 3d", c"5s 5d"]),
                            nucleus, verbosity=4)
            end
        end
    end

    @testset "Hartree–Fock equations" begin
        @testset "Hydrogen" begin
            nucleus = pc"H"
            R = RadialDifferences(N, ρ)
            atom = Atom(R, csfs(c"1s"), nucleus)

            # Set up energy expression manually.
            eng = DiagonalIntegral(o"1s")
            equations = map(Atoms.hf_equations(first(atom.configurations), eng)) do (orb,equation)
                Atoms.HFEquation(atom, equation, orb)
            end

            fock = Fock(atom, equations)

            @test length(fock.equations) == 1
            @test isapprox(SCF.energy(fock.equations[1]), -0.5, atol=2e-5)
        end

        @testset "Helium" begin
            nucleus = pc"He"
            R = RadialDifferences(N, ρ, charge(nucleus))
            atom = Atom(R, spin_configurations(c"1s2"), nucleus, verbosity=Inf)

            fock = Fock(atom; verbosity=Inf)
            scf!(fock, verbosity=Inf, num_printouts=typemax(Int))
            display(fock)
            println()

            @test all(isapprox.(SCF.energy.(fock.equations), -24.587387/27.211, rtol=0.1))
        end
    end
end
