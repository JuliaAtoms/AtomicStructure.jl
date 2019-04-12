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
            # atom = Atom(R, csfs(c"1s"), nucleus)
            atom = Atom(R, [spin_configurations(c"1s")[1]], nucleus)

            # # Set up energy expression manually.
            # eng = DiagonalIntegral(o"1s")
            # equations = map(Atoms.hf_equations(first(atom.configurations), eng)) do (orb,equation)
            #     Atoms.HFEquation(atom, equation, orb)
            # end

            # fock = Fock(atom, equations)

            fock = Fock(atom)

            @test length(fock.equations) == 1
            @test isapprox(energy(first(fock.equations)), -0.5, atol=2e-5)
        end

        @testset "Helium" begin
            nucleus = pc"He"
            R = RadialDifferences(N, ρ, charge(nucleus))
            atom = Atom(R, spin_configurations(c"1s2"), nucleus, verbosity=Inf)

            @testset "Arnoldi" begin
                fock = Fock(atom; verbosity=Inf)
                scf!(fock, verbosity=Inf, num_printouts=typemax(Int))
                display(fock)
                println()

                @test all(isapprox.(energy.(fock.equations), -24.587387/27.211, rtol=0.1))
            end

            @testset "Arnoldi shift-and-invert" begin
                fock = Fock(atom; verbosity=Inf)
                scf!(fock, verbosity=Inf, num_printouts=typemax(Int),
                     method=:arnoldi_shift_invert)
                display(fock)
                println()

                @test all(isapprox.(energy.(fock.equations), -24.587387/27.211, rtol=0.1))
            end
        end

        @testset "Beryllium" begin
            nucleus = pc"Be"
            Z = charge(nucleus)
            rₘₐₓ = 50
            ρ = 0.15
            N = ceil(Int, rₘₐₓ/ρ + 1/2)
            R = RadialDifferences(N, ρ, charge(nucleus))
            atom = Atom(R, spin_configurations(c"1s2 2s2"),
                        nucleus, verbosity=Inf)

            fock = Fock(atom)

            test_hydrogenic_slater_integrals(fock) do i
                # The energy expressions for the 1s₀α/β and 2s₂α/β
                # spin-orbitals, respectively, should be identical.
                if i ∈ 1:2
                    [(:onebody, -8, 0.07),
                     # A 1s orbital is screened by another 1s orbital and 2 2s orbitals
                     (:direct, Z*(exact_F⁰s[(o"1s",o"1s",0)] + 2exact_F⁰s[(o"2s",o"1s",0)]), 0.2),
                     (:exchange, -Z*(exact_Gᵏs[(o"2s",o"1s",0)]), 0.2)]
                elseif i ∈ 3:4
                    [(:onebody, -2, 0.01),
                     # A 2s orbital is screened by 2 1s orbitals and another 2s orbital
                     (:direct, Z*(2exact_F⁰s[(o"2s",o"1s",0)] + exact_F⁰s[(o"2s",o"2s",0)]), 0.2),
                     (:exchange, -Z*(exact_Gᵏs[(o"2s",o"1s",0)]), 0.2)]
                end
            end

            scf!(fock, ω=0.1, verbosity=2, num_printouts=typemax(Int))
            display(fock)
            println()

            @test isapprox(27.211energy(fock.equations.equations[3]), -9.3227, rtol=0.04)
            @test isapprox(27.211energy(fock.equations.equations[4]), -9.3227, rtol=0.04)
        end
    end
end