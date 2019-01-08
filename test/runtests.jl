using Atoms
using AtomicLevels
using FiniteDifferencesQuasi
using Test

@testset "Hydrogenic orbitals" begin
    rₘₐₓ = 300
    ρ = 0.25
    N = ceil(Int, rₘₐₓ/ρ + 1/2)

    @testset "Hydrogen" begin
        nucleus = pc"H"
        R=RadialDifferences(N, ρ)
        for method in [:arnoldi_shift_invert, :arnoldi, :eigen]
            atom = Atom(R, csfs([c"1s", c"2s", c"2p", c"3s", c"3p", c"3d"]),
                        nucleus, method=method, verbosity=4)
        end
    end
    @testset "Helium" begin
        nucleus = pc"He"
        R=RadialDifferences(N, ρ, charge(nucleus))
        atom = Atom(R, csfs([c"1s2", c"1s 2s", c"1s 2p", c"3s 3p", c"4s 3d", c"5s 5d"]),
                    nucleus, verbosity=4)
    end
end
