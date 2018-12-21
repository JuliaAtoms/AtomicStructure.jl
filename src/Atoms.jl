module Atoms

using AtomicLevels
import AtomicLevels: AbstractOrbital, HalfInteger

using ContinuumArrays
import ContinuumArrays.QuasiArrays: AbstractQuasiMatrix, MulQuasiArray
using LazyArrays

using LinearAlgebra

include("two_components.jl")

function unique_orbitals(csfs::Vector{C}) where {C<:CSF}
    map(csfs) do csf
        csf.config.orbitals
    end |> unique |> sort
end

# A channel represents the quantum numbers necessary to represtent a
# physical state of the many-electron wavefunction. Which quantum
# numbers are necessary depends on the problem under consideration
# (which symmetries/degeneracies are present).

const Channel{O<:AbstractOrbital,IT,TT} = Union{CSF{O,IT,TT},Level{O,IT,TT},State{O,IT,TT}}
const NonRelativisticChannel = Channel{Orbital,IntermediateTerm,Term}
const RelativisticChannel = Channel{RelativisticOrbital,HalfInteger,HalfInteger}

# A radial orbital is represented using basis coupled with a vector of
# expansion coefficients with respect to that basis; the basis
# implemented as an `AbstractQuasimatrix`. A collection of radial
# orbitals is instead represented using a matrix of such expansion
# coefficients, where each matrix column corresponds to a single
# radial orbital.

const RadialOrbital{T,B<:AbstractQuasiMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractVector}}}
const RadialOrbitals{T,B<:AbstractQuasiMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractMatrix}}}

# An atom constitutes a set of single-electron radial orbitals, CSFs
# which are comprised of anti-symmetrized combinations of such
# orbitals, and channels which are linear combinations of CSFs. The
# expansion coefficients of said channels are stored in a matrix,
# where each row corresponds to a channel and each column to CSF. In
# case the channels are CSFs, this matrix is the identity matrix.

mutable struct Atom{T,ΦT<:Union{T,TwoComponent{T}},B<:AbstractQuasiMatrix,TC<:CSF,C<:Channel,CM<:AbstractMatrix{T}}
    radial_orbitals::RadialOrbitals{ΦT,B}
    csfs::Vector{TC}
    channels::Vector{C}
    coeffs::CM
end

Atom(radial_orbitals::RadialOrbitals{ΦT,B},csfs::Vector{TC}) where {T<:Number,ΦT<:Union{T,<:TwoComponent{T}},B,TC} =
    Atom{T,ΦT,B,TC,TC,Diagonal{T}}(radial_orbitals,csfs,csfs,Diagonal(ones(T,length(csfs))))

function Atom(::Type{ΦT}, R::B, csfs::Vector{TC}) where {T<:Number,ΦT<:Union{T,TwoComponent{T}},B<:AbstractQuasiMatrix{T},TC}
    orbs = unique_orbitals(csfs)
    Φ = Matrix{ΦT}(undef, size(R,2), length(orbs))
    RΦ = MulQuasiArray{ΦT,2}(Mul(R,Φ))
    Atom(RΦ, csfs)
end

const DiracAtom{T,B,TC<:RelativisticCSF,C<:RelativisticChannel,CM} = Atom{T,TwoComponent{T},B,TC,C,CM}

Atom(R::B, csfs::Vector{TC}) where {T,B<:AbstractQuasiMatrix{T},TC} = Atom(T, R, csfs)
DiracAtom(R::B, csfs::Vector{<:TC}) where {T,B<:AbstractQuasiMatrix{T},TC<:RelativisticCSF} =
    Atom(TwoComponent{T}, R, csfs)

function Base.show(io::IO, ::MIME"text/plain", atom::Atom{T,ΦT,B,TC,C,M}) where {T,ΦT,B,TC,C,M}
    write(io, "Atom{$(ΦT),$(B)} with ")
    TC != C && write(io, "$(length(atom.channels)) $(TC)s expanded over ")
    write(io, "$(length(atom.csfs)) CSF")
    if length(atom.csfs) > 1
        write(io, "s:\n")
        show(io, "text/plain", atom.csfs)
    else
        write(io, ": ")
        show(io, atom.csfs[1])
    end
end

# The idea behind the `convert` functions below is the following: it
# is possible to convert along the chain CSF -> Level -> State (but
# not the other direction), and at each step, the expansion
# coefficients with respect to the CSFs are repeated for all channels
# C′ that are generate from an original channel C (e.g. the CSF
# `1s(₁²S|²S) 2p(₁²Pᵒ|³Pᵒ)-` generates levels with J=0,1,2).

span(::Type{Level}, csf::CSF) = levels(csf)
span(::Type{State}, level::Level) = states(level)

for (C,C′) in [(CSF,Level),(Level,State)]
    @eval begin
        function Base.convert(::Type{Atom{T,B,C′}}, atom::Atom{T,B,C,CM}) where {T,B,C<:$C,C′<:$C′,CM}
            channels′ = Vector{C′}()
            coeffs′ = Vector{typeof(atom.coeffs[1,:])}()
            for (i,c) in enumerate(atom.channels)
                c′ = span(C′, c)
                append!(channels′, c′)
                append!(coeffs′, repeat(atom.coeffs[i,:], length(c′)))
            end
            Atom(atom.radial_orbitals,atom.csfs,channels′,vcat(coeffs′...))
        end
    end
end

Base.convert(::Type{Atom{T,B,TC,C′}}, atom::Atom{T,B,TC,TC}) where {T,B,O,IT,TT,TC<:CSF{O,IT,TT},C′<:State{O,IT,TT}} =
    convert(Atom{T,B,TC,C′}, convert(Atom{T,B,Level{O,IT,TT}}, atom))

export Atom, DiracAtom

end # module
