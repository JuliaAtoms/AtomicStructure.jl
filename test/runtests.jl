using Atoms
using AtomicLevels
using AtomicPotentials
using AngularMomentumAlgebra
using SCF
# Needed as /an/ implementation of an AbstractQuasiMatrix, not the
# only possibility:
using FiniteDifferencesQuasi
using Test
using PrettyTables

include("exact_slater.jl")

function test_hydrogenic_slater_integrals!(data, tests,
                                           eq, term::Symbol, expected;
                                           print_res=false, kwargs...)
    eng = SCF.energy(eq, term)
    expectedf = convert(Float64, expected)
    δ = eng-expectedf
    test = :(isapprox($eng, $expected; $kwargs...))
    pass = @eval $test
    push!(data, [eq.orbital Dict(:onebody => "⟨h⟩", :direct => "⟨J⟩", :exchange => "⟨K⟩")[term] expected expectedf eng δ δ/abs(1e-10+abs(expectedf)) pass])
    push!(tests, test)
end

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

        @testset "Beryllium" begin
            nucleus = pc"Be"
            Z = charge(nucleus)
            R = RadialDifferences(N, ρ, charge(nucleus))
            atom = Atom(R, spin_configurations(c"1s2 2s2"),
                        nucleus, verbosity=Inf)

            fock = Fock(atom)
            data = []
            tests = []

            # The energy expressions for the 1s₀α/β and 2s₂α/β
            # spin-orbitals, respectively, should be identical.
            for eq in fock.equations[1:2]
                test_hydrogenic_slater_integrals!(data, tests, eq, :onebody, -8, rtol=0.07)
                # A 1s orbital is screened by another 1s orbital and 2 2s orbitals
                test_hydrogenic_slater_integrals!(data, tests, eq, :direct, Z*(exact_F⁰s[(o"1s",o"1s",0)] + 2exact_F⁰s[(o"2s",o"1s",0)]), rtol=0.2)
                test_hydrogenic_slater_integrals!(data, tests, eq, :exchange, -Z*(exact_Gᵏs[(o"2s",o"1s",0)]), rtol=0.2)
            end
            for eq in fock.equations[3:4]
                test_hydrogenic_slater_integrals!(data, tests, eq, :onebody, -2, rtol=0.01)
                # A 2s orbital is screened by 2 1s orbitals and another 2s orbital
                test_hydrogenic_slater_integrals!(data, tests, eq, :direct, Z*(2exact_F⁰s[(o"2s",o"1s",0)] + exact_F⁰s[(o"2s",o"2s",0)]), rtol=0.2)
                test_hydrogenic_slater_integrals!(data, tests, eq, :exchange, -Z*(exact_Gᵏs[(o"2s",o"1s",0)]), rtol=0.2)
            end
            println("Orbital energies pre-optimization:")
            # The type signature of Highlighter has changed on master of PrettyTables
            pass = Highlighter(f = (data,i,j)-> j == 8 && data[i,j], bold = true, color = :green)
            fail = Highlighter(f = (data,i,j)-> j == 8 && !data[i,j], bold = true, color = :red)
            pretty_table(vcat(data...), ["Orbital", "Term", "Expected", "Expected", "Actual", "Error", "Relative error", "Pass"],
                         highlighters=(pass,fail))
            for test in tests
                @eval @test $test
            end
            # scf!(fock, verbosity=Inf, num_printouts=typemax(Int))
            # display(fock)
            println()
        end
    end
end
