module Atoms

using AtomicLevels
import AtomicLevels: AbstractOrbital, HalfInteger
using AngularMomentumAlgebra

using ContinuumArrays
import ContinuumArrays.QuasiArrays: AbstractQuasiMatrix, MulQuasiArray, QuasiAdjoint
using LazyArrays
import LazyArrays: ⋆, materialize, materialize!, MulAdd

using LinearAlgebra
using SparseArrays

using ArnoldiMethod
using SCF
import SCF: norm_rot!, update!, KrylovWrapper, print_block
using CoulombIntegrals

using AtomicPotentials

using Formatting
using UnicodeFun

function unique_orbitals(configurations::Vector{C}) where {O,C<:Configuration{O}}
    map(configurations) do config
        config.orbitals
    end |> o -> Vector{O}(vcat(o...)) |> unique |> sort
end

getn(orb::Orbital) = orb.n
getn(orb::SpinOrbital) = orb.orb.n

getℓ(orb::Orbital) = orb.ℓ
getℓ(orb::SpinOrbital) = orb.orb.ℓ

include("radial_orbitals.jl")
include("atom_types.jl")
include("one_body.jl")
include("utils.jl")
include("hydrogenic.jl")
include("projectors.jl")
include("lazy_matmuls.jl")
include("hf_operators.jl")
include("scf_equations.jl")

end # module
