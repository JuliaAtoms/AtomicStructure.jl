# A many-electron wave function configuration can either be given as a
# CSF (summed over spins) or a configuration of spin-orbitals (where
# all quantum numbers are specified).
const ManyElectronWavefunction{O<:AbstractOrbital} = Union{CSF{O},Configuration{<:SpinOrbital}}

# An atom constitutes a set of single-electron radial orbitals,
# ManyElectronWavefunction:s which are comprised of anti-symmetrized
# combinations of such orbitals, and channels which are linear
# combinations of ManyElectronWavefunction:s. The expansion
# coefficients of said channels are stored in a matrix, where each row
# corresponds to a channel and each column to
# ManyElectronWavefunction. In case the channels are
# ManyElectronWavefunction:s, this matrix is the identity matrix.
#
# The potential can be used to model either the nucleus by itself (a
# point charge or a nucleus of finite extent) or the core orbitals
# (i.e. a pseudo-potential).

mutable struct Atom{T,B<:AbstractQuasiMatrix,O<:AbstractOrbital,TC<:ManyElectronWavefunction,C,P<:AbstractPotential}
    radial_orbitals::RadialOrbitals{T,B}
    orbitals::Vector{O}
    configurations::Vector{TC}
    mix_coeffs::Vector{C} # Mixing coefficients for multi-configurational atoms
    potential::P
end

get_config(config::Configuration) = config
get_config(csf::CSF) = csf.config

Atom(radial_orbitals::RadialOrbitals{T,B}, orbitals::Vector{O},
     configurations::Vector{<:TC}, potential::P, ::Type{C}) where {T<:Number,B,O,TC<:ManyElectronWavefunction,C,P} =
         Atom{T,B,O,TC,C,P}(radial_orbitals, orbitals,
                            configurations, vcat(one(T), zeros(C, length(configurations)-1)),
                            potential)

function Atom(::UndefInitializer, ::Type{T}, R::B, configurations::Vector{TC}, potential::P,
              ::Type{C}) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
    isempty(configurations) &&
        throw(ArgumentError("At least one configuration required to create an atom"))
    all(isequal(num_electrons(first(configurations))),
        map(config -> num_electrons(config), configurations)) ||
            throw(ArgumentError("All configurations need to have the same amount of electrons"))
    orbs = unique_orbitals(get_config.(configurations))
    Φ = Matrix{T}(undef, size(R,2), length(orbs))
    RΦ = MulQuasiArray{T,2}(Mul(R,Φ))
    Atom(RΦ, orbs, configurations, potential, C)
end

function Atom(init::Symbol, ::Type{T}, R::B, configurations::Vector{TC},
              potential::P, ::Type{C}; kwargs...) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
    atom = Atom(undef, T, R, configurations, potential, C)
    if init == :hydrogenic
        hydrogenic!(atom; kwargs...)
    else
        throw(ArgumentError("Unknown orbital initialization mode $(init)"))
    end
    atom
end

const DiracAtom{T,B,TC<:RelativisticCSF,C,P} = Atom{T,B,TC,C,P}

Atom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C};
     kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,C,P<:AbstractPotential} =
         Atom(init, T, R, configurations, potential, C; kwargs...)

DiracAtom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C};
          kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC<:RelativisticCSF,C,P<:AbstractPotential} =
    Atom(init, R, configurations, potential, C; kwargs...)

Atom(R::B, configurations::Vector{<:TC}, potential::P, ::Type{C}=eltype(R);
     kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,C,P<:AbstractPotential} =
    Atom(:hydrogenic, R, configurations, potential, C; kwargs...)

AtomicLevels.num_electrons(atom::Atom) =
    num_electrons(first(atom.configurations))

function Base.show(io::IO, atom::Atom{T,B,O,TC,C,P}) where {T,B,O,TC,C,P}
    Z = charge(atom.potential)
    N = num_electrons(atom)
    write(io, "Atom{$(T),$(B)}($(atom.potential); $(N) e⁻ ⇒ Q = $(Z-N)) with ")
    write(io, "$(length(atom.configurations)) $(TC)")
    length(atom.configurations) > 1 && write(io, "s")
end

function Base.show(io::IO, ::MIME"text/plain", atom::Atom{T,B,O,TC,C,P}) where {T,B,O,TC,C,P}
    show(io, atom)
    if length(atom.configurations) > 1
        write(io, ":\n")
        show(io, "text/plain", atom.configurations)
    else
        write(io, ": ")
        show(io, atom.configurations[1])
    end
end

radial_basis(atom::Atom) = first(atom.radial_orbitals.mul.factors)

function orbital_index(atom::Atom{T,B,O}, orb::O) where {T,B,O}
    i = findfirst(isequal(orb),atom.orbitals)
    i === nothing && throw(BoundsError("$(orb) not present among $(atom.orbitals)"))
    i
end

Base.getindex(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*atom.radial_orbitals.mul.factors[2][:,j]

Base.getindex(atom::Atom{T,B,O}, orb::O) where {T,B,O} =
    atom[orbital_index(atom, orb)]

Base.view(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*view(atom.radial_orbitals.mul.factors[2], :, j)

Base.view(atom::Atom{T,B,O}, orb::O) where {T,B,O} =
    view(atom, orbital_index(atom, orb))

"""
    norm(atom[, p=2; configuration=1])

This calculates the _amplitude_ norm of the `atom`, i.e. ᵖ√N where N
is the number electrons. By default, it uses the first `configuration`
of the `atom` to weight the individual orbital norms.
"""
function LinearAlgebra.norm(atom::Atom{T}, p::Real=2; configuration::Int=1) where T
    RT = real(T)
    n = zero(RT)
    for (orb,occ,state) in get_config(atom.configurations[configuration])
        n += occ*norm(view(atom, orb), p)^p
    end
    n^(one(RT)/p) # Unsure why you'd ever want anything but the 2-norm, but hey
end

LinearAlgebra.normalize!(fock::Fock{A}, v::V) where {A<:Atom,V<:AbstractVector} =
    normalize!(radial_basis(fock.quantum_system)*v)

function SCF.norm_rot!(ro::RO) where {RO<:RadialOrbital}
    normalize!(ro)
    SCF.rotate_max_lobe!(ro.mul.factors[2])
    ro
end

export Atom, DiracAtom, radial_basis
