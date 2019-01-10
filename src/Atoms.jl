module Atoms

using AtomicLevels
import AtomicLevels: AbstractOrbital, HalfInteger
using AngularMomentumAlgebra

using ContinuumArrays
import ContinuumArrays.QuasiArrays: AbstractQuasiMatrix, MulQuasiArray, QuasiAdjoint
using LazyArrays

using LinearAlgebra
using SparseArrays

using ArnoldiMethod
using SCF

using Formatting

function unique_orbitals(csfs::Vector{C}) where {C<:CSF}
    map(csfs) do csf
        csf.config.orbitals
    end |> o -> vcat(o...) |> unique |> sort
end

include("two_components.jl")
include("table_of_elements.jl")
include("potentials.jl")
include("atom_types.jl")
include("one_body.jl")
include("utils.jl")
include("hydrogenic.jl")
include("scf_equations.jl")

end # module