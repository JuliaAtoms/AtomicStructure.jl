exact_energies = Dict(pc"He" => [-2.8616800, [o"1s" => -0.91795555]],
                      pc"Be" => [-14.573023, [o"1s" => -4.7326698, o"2s" => -0.3092695]],

                      pc"Ne" => [-128.54710, [o"1s" => -32.7724455, o"2s" => -1.93039095, o"2p" => -0.85040965]],
                      PseudoPotentials.NeonHF => [-34.709465, [o"2s" => -1.93039095, o"2p" => -0.85040965]],
                      PseudoPotentials.NeonWB => [-34.708785, [o"2s" => -1.93039095, o"2p" => -0.85040965]],

                      pc"Mg" => [-199.61463, [o"1s" => -49.0317255, o"2s" => -3.767718, o"2p" => -2.2822236, o"3s" => -0.25305275]],

                      pc"Ar" => [-526.81751, [o"1s" => -237.22070/2, o"2s" => -24.644306/2, o"2p" => -19.142932/2,
                                             o"3s" => -2.5547063/2, o"3p" => -1.1820348/2]],
                      PseudoPotentials.ArgonHF => [-20.849867, [o"3s" => -2.5547063/2, o"3p" => -1.1820348/2]],
                      PseudoPotentials.ArgonWB => [-20.884584, [o"3s" => -2.5547063/2, o"3p" => -1.1820348/2]],

                      pc"Ca" => [-676.75818, [o"1s" => -298.72744/2, o"2s" => -33.645481/2, o"2p" => -27.258531/2,
                                          o"3s" => -4.4907488/2, o"3p" => -2.6814114/2, o"4s" => -0.3910594/2]],
                      pc"Zn" => [-1777.8481, [o"1s" => -706.60909/2, o"2s" => -88.723452/2, o"2p" => -77.849691/2,
                                             o"3s" => -11.275642/2, o"3p" => -7.6787581/2, o"3d" => -1.5650843/2,
                                             o"4s" => -0.5850141/2]],
                      pc"Kr" => [-2752.0550, [o"1s" => -1040.3309/2, o"2s" => -139.80617/2, o"2p" => -126.01957/2,
                                             o"3s" => -21.698934/2, o"3p" => -16.663004/2, o"3d" => -7.6504697/2,
                                             o"4s" => -2.3058703/2, o"4p" => -1.0483734/2]],
                      pc"Xe" => [-7232.1384, [o"1s" => -2448.7956/2, o"2s" => -378.68024/2, o"2p" => -355.56490/2,
                                             o"3s" => -80.351324/2, o"3p" => -70.443321/2, o"3d" => -52.237736/2,
                                             o"4s" => -15.712602/2, o"4p" => -12.016674, o"4d" => -5.5557597,
                                             o"5s" => -1.8888279/2, o"5p" => -0.9145793/2]],

                      PseudoPotentials.XenonHF => [-14.989100, [o"5s" => -1.8888279/2, o"5p" => -0.9145793/2]],
                      PseudoPotentials.XenonWB => [-15.277055, [o"5s" => -1.8888279/2, o"5p" => -0.9145793/2]])

function energy_errors(fock, exact_energies, Δ, δ)
    atom = fock.quantum_system

    Eexact,orbExact = exact_energies
    exact_energies = vcat(Eexact, [repeat([E], degeneracy(o)) for (o,E) in orbExact]...)

    H = zeros(1,1)
    Etot = SCF.energy_matrix!(H, fock)[1,1]
    energies = vcat(Etot, SCF.energy.(fock.equations))
    errors = energies - exact_energies

    labels = vcat("Total", string.(atom.orbitals))
    pretty_table([labels energies 27.211energies exact_energies errors 27.211errors errors./abs.(exact_energies)],
                 ["", "Energy [Ha]", "Energy [eV]", "HF limit [Ha]", "Error [Ha]", "Error [eV]", "Relative error"])

    @test abs(errors[1]) < Δ
    @test all(abs.(errors[2:end]) .< δ)
end

function atom_calc(nucleus::AbstractPotential, grid_type::Symbol, rₘₐₓ, ρ,
                   Δ, δ; max_iter=200, kwargs...)
    # If we're using a pseudopotential, we don't really need higher
    # order of the FEDVR basis functions in the first finite element,
    # since the orbitals almost vanish there.
    R,r = get_atom_grid(grid_type, rₘₐₓ, ρ, nucleus,
                        amend_order=nucleus isa PointCharge)

    atom = Atom(R, [spin_configurations(ground_state(nucleus))[1]],
                nucleus, eltype(R))

    fock = Fock(atom)

    scf!(fock; max_iter=max_iter, num_printouts=typemax(Int), verbosity=2, kwargs...)
    energy_errors(fock, exact_energies[nucleus], Δ, δ)
end

@testset "Calculation accuracy" begin
    @testset "Helium" begin
        @testset "$(grid_type)" for (grid_type,Δ,δ) in [(:fedvr,2e-8,1e-9),
                                                        (:fd,4e-3,4e-3)]
            atom_calc(pc"He", grid_type, 10.0, 0.1, Δ, δ, ω=0.9)
        end
    end

    @testset "Beryllium" begin
        @testset "$(grid_type)" for (grid_type,Δ,δ) in [(:fedvr,2e-8,0.02),
                                                        (:fd,0.06,0.025)]
            atom_calc(pc"Be", grid_type, 15.0, 0.1, Δ, δ, ω=0.9)
        end
    end

    @testset "Neon" begin
        @testset "$nucleus" for (nucleus,Δ,δ) in [(pc"Ne",0.005,0.05),
                                                  (PseudoPotentials.NeonHF,0.2,0.05)]
            atom_calc(nucleus, :fedvr, 10, 0.2, Δ, δ,
                      ω=0.9, ωmax=1.0-1e-5)
        end
    end

    @testset "Xenon" begin
        @testset "$nucleus" for (nucleus,Δ,δ) in [# (pc"Xe",0.005,0.05),
                                                  (PseudoPotentials.XenonHF,0.005,0.005)]
            atom_calc(nucleus, :fedvr, 10, 0.2, Δ, δ,
                      ω=0.9, ωmax=1.0-1e-5)
        end
    end
end
