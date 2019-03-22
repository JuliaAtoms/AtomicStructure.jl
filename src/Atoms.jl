module Atoms

using AtomicLevels
import AtomicLevels: AbstractOrbital, HalfInteger
using AngularMomentumAlgebra
using EnergyExpressions
import EnergyExpressions: NBodyMatrixElement, OrbitalMatrixElement,
    orbital_equation, MCEquationSystem, QuantumOperator

using ContinuumArrays
import ContinuumArrays: Basis
import ContinuumArrays.QuasiArrays: AbstractQuasiMatrix, MulQuasiArray, QuasiAdjoint
using LazyArrays
import LazyArrays: ⋆, materialize, materialize!, MulAdd
using FillArrays
using BandedMatrices

using FiniteDifferencesQuasi

using LinearAlgebra
using SparseArrays

using ArnoldiMethod
using SCF
import SCF: norm_rot!, update!, KrylovWrapper, print_block
using CoulombIntegrals
using MatrixFactorizations

using AtomicPotentials

using Formatting
using UnicodeFun

function unique_orbitals(configurations::Vector{C}) where {O,C<:Configuration{O}}
    map(configurations) do config
        config.orbitals
    end |> o -> Vector{O}(vcat(o...)) |> unique |> sort
end

getspatialorb(orb::Orbital) = orb
getspatialorb(orb::SpinOrbital) = orb.orb

getn(orb::AbstractOrbital) = getspatialorb(orb).n
getℓ(orb::AbstractOrbital) = getspatialorb(orb).ℓ

include("restricted_bases.jl")
include("radial_orbitals.jl")
include("atom_types.jl")
include("one_body.jl")
include("utils.jl")
include("hydrogenic.jl")
include("projectors.jl")
include("lazy_matmuls.jl")
include("orbital_integrals.jl")
include("orbital_hamiltonian.jl")
include("orbital_equation.jl")
include("common_integrals.jl")
include("observables.jl")
include("equations.jl")

end # module
