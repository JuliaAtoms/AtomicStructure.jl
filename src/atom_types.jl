"""
    ManyElectronWavefunction

A many-electron wave function configuration can either be given as a
CSF (summed over spins) or a configuration of spin-orbitals (where all
quantum numbers are specified).
"""
const ManyElectronWavefunction{O<:AbstractOrbital} = Union{CSF{O},Configuration{<:SpinOrbital}}

"""
    Atom(radial_orbitals, orbitals, configurations, mix_coeffs, potential)

An atom constitutes a set of single-electron `orbitals` with
associated `radial_orbitals`, `configurations` which are
[`ManyElectronWavefunction`](@ref):s, comprising of anti-symmetrized
combinations of such orbitals. The expansion coefficients `mix_coeffs`
determine the linear combination of the `configurations` for
multi-configurational atoms.

The `potential` can be used to model either the nucleus by itself (a
point charge or a nucleus of finite extent) or the core orbitals
(i.e. a pseudo-potential).
"""
mutable struct Atom{T,B<:AbstractQuasiMatrix,O<:AbstractOrbital,TC<:ManyElectronWavefunction,C,P<:AbstractPotential} <: AbstractQuantumSystem
    radial_orbitals::RadialOrbitals{T,B}
    orbitals::Vector{O}
    configurations::Vector{TC}
    mix_coeffs::Vector{C} # Mixing coefficients for multi-configurational atoms
    potential::P
end

get_config(config::Configuration) = config
get_config(csf::CSF) = csf.config

"""
    outsidecoremodel(configuration::Configuration, potential::P)

Return the part of the electronic `configuration` that is not part of
the the configuration modelled by the `potential`. For a point charge,
this is the same as the `configuration` itself, but for
pseudopotential, typically only the outer shells remain.
"""
outsidecoremodel(configuration::Configuration, potential::P) where P =
    filter((orb,occ,state) -> getspatialorb(orb) ∉ core(ground_state(potential)), configuration)
outsidecoremodel(configuration::Configuration, ::PointCharge) where P =
    configuration

"""
    Atom(radial_orbitals, orbitals, configurations, potential, ::Type{C})

Create an [`Atom`](@ref) from the lists of `radial_orbitals`
associated with `orbitals` and electronic `configurations`, with a
nucleus modelled by `potential`. `C` determines the `eltype` of the
mixing coefficients, which are initialized to `[1,0,0,...]`.
"""
Atom(radial_orbitals::RadialOrbitals{T,B}, orbitals::Vector{O},
     configurations::Vector{<:TC}, potential::P, ::Type{C}) where {T<:Number,B,O,TC<:ManyElectronWavefunction,C,P} =
         Atom{T,B,O,TC,C,P}(radial_orbitals, orbitals,
                            configurations, vcat(one(T), zeros(C, length(configurations)-1)),
                            potential)

"""
    Atom(undef, ::Type{T}, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})

Create an [`Atom`](@ref) on the space spanned by `R`, from the list of
electronic `configurations`, with a nucleus modelled by `potential`,
and leave the orbitals uninitialized. `T` determines the `eltype` of
the radial orbitals and `C` the mixing coefficients.
"""
function Atom(::UndefInitializer, ::Type{T}, R::B, configurations::Vector{TC}, potential::P,
              ::Type{C}) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
    isempty(configurations) &&
        throw(ArgumentError("At least one configuration required to create an atom"))
    all(isequal(num_electrons(first(configurations))),
        map(config -> num_electrons(config), configurations)) ||
            throw(ArgumentError("All configurations need to have the same amount of electrons"))
    all(isequal(core(get_config(first(configurations)))),
        map(config -> core(get_config(config)), configurations)) ||
            throw(ArgumentError("All configurations need to share the same core orbitals"))

    pot_cfg = core(ground_state(potential))
    for cfg in configurations
        cfg_closed_orbitals = getspatialorb.(core(get_config(cfg)).orbitals)
        pot_cfg.orbitals ⊆ cfg_closed_orbitals ||
            throw(ArgumentError("Configuration modelled by nuclear potential ($(pot_cfg)) must belong to the closed set of all configurations"))
    end

    orbs = unique_orbitals(outsidecoremodel.(get_config.(configurations), Ref(potential)))

    Φ = Matrix{T}(undef, size(R,2), length(orbs))
    RΦ = MulQuasiArray{T,2}(Mul(R,Φ))
    Atom(RΦ, orbs, configurations, potential, C)
end

"""
    Atom(init, ::Type{T}, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})

Create an [`Atom`](@ref) on the space spanned by `R`, from the list of
electronic `configurations`, with a nucleus modelled by `potential`,
and initialize the orbitals according to `init`. `T` determines the
`eltype` of the radial orbitals and `C` the mixing coefficients.
"""
function Atom(init::Symbol, ::Type{T}, R::B, configurations::Vector{TC},
              potential::P, ::Type{C}; kwargs...) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
    atom = Atom(undef, T, R, configurations, potential, C)
    if init == :hydrogenic
        hydrogenic!(atom; kwargs...)
    elseif init == :zero
        atom.radial_orbitals.mul.factors[2] .= zero(T)
    else
        throw(ArgumentError("Unknown orbital initialization mode $(init)"))
    end
    atom
end

"""
    DiracAtom

A `DiracAtom` is a specialization of [`Atom`](@ref) for the
relativistic case.
"""
const DiracAtom{T,B,TC<:RelativisticCSF,C,P} = Atom{T,B,TC,C,P}

"""
    Atom(init, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})

Create an [`Atom`](@ref) on the space spanned by `R`, from the list of
electronic `configurations`, with a nucleus modelled by `potential`,
and initialize the orbitals according to `init`. `C` determines the
`eltype` of the mixing coefficients.
"""
Atom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C};
     kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,C,P<:AbstractPotential} =
         Atom(init, T, R, configurations, potential, C; kwargs...)

DiracAtom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C};
          kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC<:RelativisticCSF,C,P<:AbstractPotential} =
    Atom(init, T, R, configurations, potential, C; kwargs...)

"""
    Atom(R::AbstractQuasiMatrix, configurations, potential[, ::Type{C}=eltype(R)])

Create an [`Atom`](@ref) on the space spanned by `R`, from the list of
electronic `configurations`, with a nucleus modelled by `potential`,
and initialize the orbitals to their hydrogenic values.
"""
Atom(R::B, configurations::Vector{<:TC}, potential::P, ::Type{C}=eltype(R);
     kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,C,P<:AbstractPotential} =
         Atom(:hydrogenic, R, configurations, potential, C; kwargs...)

"""
    Atom(other_atom::Atom, configurations)

Create a new atom using the same basis and nuclear potential as
`other_atom`, but with a different set of `configurations`. The
orbitals of `other_atom` are copied over as starting guess.
"""
function Atom(other_atom::Atom{T,B,O,TC,C,P},
              configurations::Vector{<:TC}; kwargs...) where {T,B,O,TC,C,P}
    core(first(configurations)) == core(first(other_atom.configurations)) ||
        throw(ArgumentError("Core orbitals must be the same"))
    R = radial_basis(other_atom)
    # It is slightly wasteful to first initialize all orbitals
    # hydrogenically, and the overwrite with orbitals from
    # other_atom. Ideally, missing orbitals would be estimated by
    # diagonalizing e.g. their energy expressions without the exchange
    # terms.
    atom = Atom(R, configurations, other_atom.potential, C; kwargs...)
    copyto!(atom, other_atom)
    atom
end

"""
    copyto!(dst::Atom, src::Atom)

Copy the radial orbitals of `src` to the corresponding ones in
`dst`. Orbitals of `src` missing from `dst` are skipped. It is assumed
that the underlying bases are compatible, i.e. that the radial
coordinate of `dst` encompasses that of `src` and grid spacing
&c. agree.
"""
function Base.copyto!(dst::Atom, src::Atom)
    Rsrc = radial_basis(src)
    Rdst = radial_basis(dst)
    axes(Rsrc,1) ⊆ axes(Rdst,1) &&
        axes(Rsrc,2) ⊆ axes(Rdst,2) ||
        throw(ArgumentError("$(Rsrc) not a subset of $(Rdst)"))

    for o in src.orbitals
        o ∈ dst.orbitals || continue
        copyto!(view(dst, o).mul.factors[2],
                view(src, o).mul.factors[2])
    end
end

"""
    num_electrons(atom)

Return number of electrons in `atom`.
"""
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

"""
    getindex(atom, j)

Returns a copy of the `j`:th radial orbital.
"""
Base.getindex(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*atom.radial_orbitals.mul.factors[2][:,j]

"""
    getindex(atom, orb)

Returns a copy of the radial orbital corresponding to `orb`.
"""
Base.getindex(atom::Atom{T,B,O}, orb::O) where {T,B,O} =
    atom[orbital_index(atom, orb)]

"""
    view(atom, j)

Returns a `view` of the `j`:th radial orbital.
"""
Base.view(atom::Atom, j::I) where {I<:Integer} =
    radial_basis(atom)*view(atom.radial_orbitals.mul.factors[2], :, j)

"""
    view(atom, orb)

Returns a `view` of the radial orbital corresponding to `orb`.
"""
Base.view(atom::Atom{T,B,O}, orb::O) where {T,B,O} =
    view(atom, orbital_index(atom, orb))

"""
    SCF.coefficients(atom)

Returns a `view` of the mixing coefficients.
"""
SCF.coefficients(atom::A) where {A<:Atom} =
    view(atom.mix_coeffs, :)

"""
    SCF.orbitals(atom)

Returns a `view` of the radial orbital coefficients (NB, it does _not_
return the `MulQuasiMatrix`, but the actual underlying expansion
coefficients, since `SCF` operates on them in the self-consistent
iteration).
"""
SCF.orbitals(atom::A) where {A<:Atom} =
    view(atom.radial_orbitals.mul.factors[2], :, :)

"""
    norm(atom[, p=2; configuration=1])

This calculates the _amplitude_ norm of the `atom`, i.e. ᵖ√N where N
is the number electrons. By default, it uses the first `configuration`
of the `atom` to weight the individual orbital norms.
"""
function LinearAlgebra.norm(atom::Atom{T}, p::Real=2; configuration::Int=1) where T
    RT = real(T)
    n = zero(RT)
    for (orb,occ,state) in outsidecoremodel(get_config(atom.configurations[configuration]),
                                            atom.potential)
        n += occ*norm(view(atom, orb), p)^p
    end
    n^(one(RT)/p) # Unsure why you'd ever want anything but the 2-norm, but hey
end

LinearAlgebra.normalize!(atom::A, v::V) where {A<:Atom,V<:AbstractVector} =
    normalize!(radial_basis(atom)*v)

function SCF.norm_rot!(ro::RO) where {RO<:RadialOrbital}
    normalize!(ro)
    SCF.rotate_max_lobe!(ro.mul.factors[2])
    ro
end

export Atom, DiracAtom, radial_basis
