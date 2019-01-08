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

const RadialCoeff{T<:Number} = Union{T,TwoComponent{T}}
const RadialOrbital{T<:RadialCoeff,B<:AbstractQuasiMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractVector}}}
const RadialOrbitals{T<:RadialCoeff,B<:AbstractQuasiMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractMatrix}}}

# An atom constitutes a set of single-electron radial orbitals, CSFs
# which are comprised of anti-symmetrized combinations of such
# orbitals, and channels which are linear combinations of CSFs. The
# expansion coefficients of said channels are stored in a matrix,
# where each row corresponds to a channel and each column to CSF. In
# case the channels are CSFs, this matrix is the identity matrix.
#
# The potential can be used to model either the nucleus by itself (a
# point charge or a nucleus of finite extent) or the core orbitals
# (i.e. a pseudo-potential).

mutable struct Atom{T,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix,O<:AbstractOrbital,TC<:CSF,C<:Channel,CM<:AbstractMatrix{T},P<:AbstractPotential}
    radial_orbitals::RadialOrbitals{ΦT,B}
    orbitals::Vector{O}
    csfs::Vector{TC}
    csf_orbital_indices::SparseMatrixCSC{Int}
    channels::Vector{C}
    coeffs::CM
    potential::P
end

"""
   get_csf_orbital_indices(orbitals, csfs)

Return matrix where each row indicates how many of each orbital the
corresponding CSF is comprised of.
"""
function get_csf_orbital_indices(orbitals::Vector{O}, csfs::Vector{TC}) where {O<:AbstractOrbital,TC<:CSF}
    indices = spzeros(Int, length(csfs), length(orbitals))
    for (i,csf) in enumerate(csfs)
        for (orb,occ,state) in csf.config
            indices[i,findfirst(isequal(orb), orbitals)] = occ
        end
    end
    indices
end

Atom(radial_orbitals::RadialOrbitals{ΦT,B}, orbitals::Vector{O},
     csfs::Vector{<:TC}, potential::P) where {T<:Number,ΦT<:RadialCoeff{T},B,O,TC<:CSF,P} =
    Atom{T,ΦT,B,O,TC,TC,Diagonal{T},P}(radial_orbitals, orbitals,
                                       csfs, get_csf_orbital_indices(orbitals, csfs),
                                       csfs, Diagonal(ones(T,length(csfs))),
                                       potential)

function Atom(::UndefInitializer, ::Type{ΦT}, R::B, csfs::Vector{TC}, potential::P) where {T<:Number,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix{T},TC,P}
    orbs = unique_orbitals(csfs)
    Φ = Matrix{ΦT}(undef, size(R,2), length(orbs))
    RΦ = MulQuasiArray{ΦT,2}(Mul(R,Φ))
    Atom(RΦ, orbs, csfs, potential)
end

function Atom(init::Symbol, ::Type{ΦT}, R::B, csfs::Vector{TC},
              potential::P; kwargs...) where {T<:Number,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix{T},TC,P}
    atom = Atom(undef, ΦT, R, csfs, potential)
    if init == :hydrogenic
        hydrogenic!(atom; kwargs...)
    else
        throw(ArgumentError("Unknown orbital initialization mode $(init)"))
    end
end

const DiracAtom{T,B,TC<:RelativisticCSF,C<:RelativisticChannel,CM,P} = Atom{T,TwoComponent{T},B,TC,C,CM,P}

Atom(init::Init, R::B, csfs::Vector{TC}, potential::P; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,P<:AbstractPotential} =
    Atom(init, T, R, csfs, potential; kwargs...)
DiracAtom(init::Init, R::B, csfs::Vector{TC}, potential::P; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC<:RelativisticCSF,P<:AbstractPotential} =
    Atom(init, TwoComponent{T}, R, csfs, potential; kwargs...)

Atom(R::B, csfs::Vector{TC}, potential::P; kwargs...) where {B<:AbstractQuasiMatrix,TC<:CSF,P<:AbstractPotential} =
    Atom(:hydrogenic, R, csfs, potential; kwargs...)

function Base.show(io::IO, atom::Atom{T,ΦT,B,O,TC,C,M,P}) where {T,ΦT,B,O,TC,C,M,P}
    write(io, "Atom{$(ΦT),$(B)}($(atom.potential)) with ")
    TC != C && write(io, "$(length(atom.channels)) $(TC)s expanded over ")
    write(io, "$(length(atom.csfs)) CSF")
    length(atom.csfs) > 1 && write(io, "s")
end

function Base.show(io::IO, ::MIME"text/plain", atom::Atom{T,ΦT,B,O,TC,C,M,P}) where {T,ΦT,B,O,TC,C,M,P}
    show(io, atom)
    if length(atom.csfs) > 1
        write(io, ":\n")
        show(io, "text/plain", atom.csfs)
    else
        write(io, ": ")
        show(io, atom.csfs[1])
    end
end

radial_basis(atom::Atom) = first(atom.radial_orbitals.mul.factors)

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
        function Base.convert(::Type{Atom{T,B,O,C′}}, atom::Atom{T,B,O,C,CM}) where {T,B,O,C<:$C,C′<:$C′,CM}
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

Base.convert(::Type{Atom{T,B,O,TC,C′}}, atom::Atom{T,B,O,TC,TC}) where {T,B,O,IT,TT,TC<:CSF{O,IT,TT},C′<:State{O,IT,TT}} =
    convert(Atom{T,B,O,TC,C′}, convert(Atom{T,B,O,Level{O,IT,TT}}, atom))

export Atom, DiracAtom, radial_basis
