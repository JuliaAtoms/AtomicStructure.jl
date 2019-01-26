# A many-electron wave function configuration can either be given as a
# CSF (summed over spins) or a configuration of spin-orbitals (where
# all quantum numbers are specified).
const ManyElectronWavefunction{O<:AbstractOrbital} = Union{CSF{O},Configuration{<:SpinOrbital}}

# A channel represents the quantum numbers necessary to represent a
# physical state of the many-electron wavefunction. Which quantum
# numbers are necessary depends on the problem under consideration
# (which symmetries/degeneracies are present).

const Channel{O<:AbstractOrbital,IT,TT} = Union{CSF{O,IT,TT},Level{O,IT,TT},State{O,IT,TT},Configuration{<:SpinOrbital}}
const NonRelativisticChannel = Channel{Orbital,IntermediateTerm,Term}
const RelativisticChannel = Channel{RelativisticOrbital,HalfInteger,HalfInteger}

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

mutable struct Atom{T,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix,O<:AbstractOrbital,TC<:ManyElectronWavefunction,C<:Channel,CM<:AbstractMatrix{T},P<:AbstractPotential}
    radial_orbitals::RadialOrbitals{ΦT,B}
    orbitals::Vector{O}
    configurations::Vector{TC}
    config_orbital_indices::SparseMatrixCSC{Int}
    channels::Vector{C}
    coeffs::CM
    potential::P
end

"""
   get_config_orbital_indices(orbitals, configurations)

Return matrix where each row indicates how many of each orbital the
corresponding configuration is comprised of.
"""
function get_config_orbital_indices(orbitals::Vector{O}, configs::Vector{TC}) where {O<:AbstractOrbital,TC<:Configuration}
    indices = spzeros(Int, length(configs), length(orbitals))
    for (i,config) in enumerate(configs)
        for (orb,occ,state) in config
            indices[i,findfirst(isequal(orb), orbitals)] = occ
        end
    end
    indices
end

get_config(config::Configuration) = config
get_config(csf::CSF) = csf.config

Atom(radial_orbitals::RadialOrbitals{ΦT,B}, orbitals::Vector{O},
     configurations::Vector{<:TC}, potential::P) where {T<:Number,ΦT<:RadialCoeff{T},B,O,TC<:ManyElectronWavefunction,P} =
    Atom{T,ΦT,B,O,TC,TC,Diagonal{T},P}(radial_orbitals, orbitals,
                                       configurations, get_config_orbital_indices(orbitals, get_config.(configurations)),
                                       configurations, Diagonal(ones(T,length(configurations))),
                                       potential)

function Atom(::UndefInitializer, ::Type{ΦT}, R::B, configurations::Vector{TC}, potential::P) where {T<:Number,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix{T},TC,P}
    isempty(configurations) &&
        throw(ArgumentError("At least one configuration required to create an atom"))
    all(isequal(num_electrons(first(configurations))),
        map(config -> num_electrons(config), configurations)) ||
            throw(ArgumentError("All configurations need to have the same amount of electrons"))
    orbs = unique_orbitals(get_config.(configurations))
    Φ = Matrix{ΦT}(undef, size(R,2), length(orbs))
    RΦ = MulQuasiArray{ΦT,2}(Mul(R,Φ))
    Atom(RΦ, orbs, configurations, potential)
end

function Atom(init::Symbol, ::Type{ΦT}, R::B, configurations::Vector{TC},
              potential::P; kwargs...) where {T<:Number,ΦT<:RadialCoeff{T},B<:AbstractQuasiMatrix{T},TC,P}
    atom = Atom(undef, ΦT, R, configurations, potential)
    if init == :hydrogenic
        hydrogenic!(atom; kwargs...)
    else
        throw(ArgumentError("Unknown orbital initialization mode $(init)"))
    end
end

const DiracAtom{T,B,TC<:RelativisticCSF,C<:RelativisticChannel,CM,P} = Atom{T,TwoComponent{T},B,TC,C,CM,P}

Atom(init::Init, R::B, configurations::Vector{TC}, potential::P; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,P<:AbstractPotential} =
    Atom(init, T, R, configurations, potential; kwargs...)
DiracAtom(init::Init, R::B, configurations::Vector{TC}, potential::P; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC<:RelativisticCSF,P<:AbstractPotential} =
    Atom(init, TwoComponent{T}, R, configurations, potential; kwargs...)

Atom(R::B, configurations::Vector{<:TC}, potential::P; kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,P<:AbstractPotential} =
    Atom(:hydrogenic, R, configurations, potential; kwargs...)

AtomicLevels.num_electrons(atom::Atom) =
    num_electrons(first(atom.configurations))

function Base.show(io::IO, atom::Atom{T,ΦT,B,O,TC,C,M,P}) where {T,ΦT,B,O,TC,C,M,P}
    Z = charge(atom.potential)
    N = num_electrons(atom)
    write(io, "Atom{$(ΦT),$(B)}($(atom.potential); $(N) e⁻ ⇒ Q = $(Z-N)) with ")
    TC != C && write(io, "$(length(atom.channels)) $(C)s expanded over ")
    write(io, "$(length(atom.configurations)) $(TC)")
    length(atom.configurations) > 1 && write(io, "s")
end

function Base.show(io::IO, ::MIME"text/plain", atom::Atom{T,ΦT,B,O,TC,C,M,P}) where {T,ΦT,B,O,TC,C,M,P}
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

function orbital_index(atom::Atom{T,ΦT,B,O}, orb::O) where {T,ΦT,B,O}
    i = findfirst(isequal(orb),atom.orbitals)
    i === nothing && throw(BoundsError("$(orb) not present among $(atom.orbitals)"))
    i
end

Base.getindex(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*atom.radial_orbitals.mul.factors[2][:,j]

Base.getindex(atom::Atom{T,ΦT,B,O}, orb::O) where {T,ΦT,B,O} =
    atom[orbital_index(atom, orb)]

Base.view(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*view(atom.radial_orbitals.mul.factors[2], :, j)

Base.view(atom::Atom{T,ΦT,B,O}, orb::O) where {T,ΦT,B,O} =
    view(atom, orbital_index(atom, orb))

"""
    norm(atom[, p=2; configuration=1])

This calculates the _amplitude_ norm of the `atom`, i.e. ᵖ√N where N
is the number electrons. By default, it uses the first `configuration`
of the `atom` to weight the individual orbital norms.
"""
function LinearAlgebra.norm(atom::Atom{T}, p::Real=2; configuration::Int=1) where T
    n = zero(T)
    for (orb,occ,state) in get_config(atom.configurations[configuration])
        n += occ*norm(view(atom, orb), p)^p
    end
    n^(one(T)/p) # Unsure why you'd ever want anything but the 2-norm, but hey
end

LinearAlgebra.normalize!(fock::Fock{A}, v::V) where {A<:Atom,V<:AbstractVector} =
    normalize!(radial_basis(fock.quantum_system)*v)

function SCF.norm_rot!(ro::RO) where {RO<:RadialOrbital}
    normalize!(ro)
    SCF.rotate_max_lobe!(ro.mul.factors[2])
    ro
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
        function Base.convert(::Type{Atom{T,B,O,C′}}, atom::Atom{T,B,O,C,CM}) where {T,B,O,C<:$C,C′<:$C′,CM}
            channels′ = Vector{C′}()
            coeffs′ = Vector{typeof(atom.coeffs[1,:])}()
            for (i,c) in enumerate(atom.channels)
                c′ = span(C′, c)
                append!(channels′, c′)
                append!(coeffs′, repeat(atom.coeffs[i,:], length(c′)))
            end
            Atom(atom.radial_orbitals,atom.configurations,channels′,vcat(coeffs′...))
        end
    end
end

Base.convert(::Type{Atom{T,B,O,TC,C′}}, atom::Atom{T,B,O,TC,TC}) where {T,B,O,IT,TT,TC<:CSF{O,IT,TT},C′<:State{O,IT,TT}} =
    convert(Atom{T,B,O,TC,C′}, convert(Atom{T,B,O,Level{O,IT,TT}}, atom))

export Atom, DiracAtom, radial_basis
