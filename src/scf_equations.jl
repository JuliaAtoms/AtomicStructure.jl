# * Hartree–Fock equation

mutable struct HFEquation{T, ΦT, #<:RadialCoeff{T},
                          B<:AbstractQuasiMatrix,
                          O<:AbstractOrbital,
                          A<:Atom{T,ΦT,B,O},
                          E,
                          M# <:OrbitalHamiltonian{T,ΦT,B,O}
                          }
    atom::A
    equation::E
    orbital::O
    ϕ::RadialOrbital{T,B}
    hamiltonian::M
end

SCF.update!(eq::HFEquation; kwargs...) = update!(eq.hamiltonian; kwargs...)
SCF.KrylovWrapper(eq::HFEquation) = KrylovWrapper(eq.hamiltonian)

const EnergyExpression = Tuple{Vector{<:OneBodyHamiltonian},Tuple{Vector{<:DirectExchangePotentials},Vector{NTuple{2,Vector{Pair{Int,Float64}}}}}}

function orbital_hamiltonian(atom::A, (one_body,(two_body,multipole_terms))::E,
                             orbital::O₁,
                             projector::Proj) where {T,ΦT, #<:RadialCoeff{T},
                                                     B<:AbstractQuasiMatrix,
                                                     O₁<:AbstractOrbital,
                                                     O₂<:AbstractOrbital,
                                                     A<:Atom{T,ΦT,B,O₂},
                                                     E<:EnergyExpression,
                                                     Proj}
    R = radial_basis(atom)

    OV = typeof(view(atom, 1))
    O = promote_type(O₁, O₂)
    PO = typeof(R*Diagonal(Vector{T}(undef, size(R,2)))*R')
    direct_potentials = Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}()
    exchange_potentials = Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}()

    count(.!iszero.(one_body)) ≤ 1 ||
        throw(ArgumentError("There can only be one one-body Hamiltonian per orbital"))

    ĥ = if isempty(one_body)
        zero(AtomicOneBodyHamiltonian(atom, orbital))
    else
        AtomicOneBodyHamiltonian(atom, one_body[1].orb)
    end

    action_orbital = if !isempty(one_body)
        one_body[1].orb
    elseif !isempty(two_body)
        two_body[1].o
    else
        orbital
    end

    for (tb,mpt) ∈ zip(two_body,multipole_terms)
        a = tb.a
        b = tb.b
        av = view(atom, a)
        bv = view(atom, b)

        for (k,c) in mpt[1]
            push!(direct_potentials,
                  HFPotential(:direct, k, a, b, av, bv) => c)
        end
        for (k,c) in mpt[2]
            push!(exchange_potentials,
                  HFPotential(:exchange, k, a, b, av, bv) => c)
        end
    end

    OrbitalSplitHamiltonian(R, ĥ, direct_potentials, exchange_potentials,
                            projector,action_orbital)
end

function HFEquation(atom::A, (one_body,(two_body,multipole_terms))::E,
                    orbital::O,
                    symmetry_orbitals::Vector{O}) where {T,ΦT, #<:RadialCoeff{T},
                                                         B<:AbstractQuasiMatrix,
                                                         O<:AbstractOrbital,
                                                         A<:Atom{T,ΦT,B,O},
                                                         E<:EnergyExpression}
    OV = typeof(view(atom, 1))
    projector = Projector(OV[view(atom, other)
                             for other in symmetry_orbitals])

    hamiltonian = orbital_hamiltonian(atom, (one_body,(two_body,multipole_terms)),
                                      orbital, projector)

    HFEquation(atom, (one_body,two_body), orbital, view(atom, orbital), hamiltonian)
end

SCF.energy(hfeq::HFEquation{E,O,M}, term::Symbol=:all) where {E,O,M} = (hfeq.ϕ' * hfeq.hamiltonian[term] * hfeq.ϕ)[1]

function Base.show(io::IO, hfeq::HFEquation)
    write(io, "Hartree–Fock equation: 0 = [𝓗  - E($(hfeq.orbital))]|$(hfeq.orbital)⟩ = $(hfeq.hamiltonian) - E($(hfeq.orbital))|$(hfeq.orbital)⟩")
    EHa = SCF.energy(hfeq)
    write(io, "\n    ⟨$(hfeq.orbital)| 𝓗 |$(hfeq.orbital)⟩ = $(EHa) Ha = $(27.211EHa) eV")
    Eh = (hfeq.ϕ' * hfeq.hamiltonian[:onebody] * hfeq.ϕ)[1]
    Ed = (hfeq.ϕ' * hfeq.hamiltonian[:direct] * hfeq.ϕ)[1]
    Ex = (hfeq.ϕ' * hfeq.hamiltonian[:exchange] * hfeq.ϕ)[1]
    write(io, " (⟨h⟩ = $(Eh) Ha, ⟨J⟩ = $(Ed) Ha, ⟨K⟩ = $(Ex) Ha)")
end

# * Orbital symmetries
"""
    find_symmetries(orbitals)

Group all orbitals according to their symmetries, e.g. ℓ for
`Orbital`s. This is used to determine which off-diagonal Lagrange
multipliers are necessary to maintain orthogonality.
"""
find_symmetries(orbitals::Vector{O}) where {O<:AbstractOrbital} =
    merge!(vcat, [Dict(symmetry(orb) => [orb])
                  for orb in orbitals]...)

# * Setup Hartree–Fock equations
#=
function hf_equations(csf::NonRelativisticCSF, eng::EnergyExpression; verbosity=0)
    pconfig = peel(csf.config)
    orbitals = pconfig.orbitals

    λs = lagrange_multipliers(orbitals)
    eng += λs

    if verbosity > 2
        println("Orbitals: $(orbitals)")
        println("Lagrange multipliers: $(λs)")
    end
    verbosity > 1 && println("E = $eng\n")

    map(pconfig) do (orb,occ,state)
        equation = diff(eng, orb, occ)
        verbosity > 2 && println("\n∂[$(orb)] E = $(equation)")
        orb => equation
    end
end

"""
    hf_equations(csf[; verbosity=0])

Derive the Hartree–Fock equations for the non-relativistic
configuration state function `csf`. All constituent orbitals are
assumed spectroscopic, i.e. they are all assigned Lagrange multipliers
to ensure their orthonormality. Additionally, orbitals of the same `ℓ`
but different `n` are assigned off-diagonal Lagrange multipliers.
"""
function hf_equations(csf::NonRelativisticCSF; verbosity=0)
    verbosity > 0 &&
        println("Deriving HF equations for $(csf)")

    # energy_expression has to be provided by an angular momentum
    # library.
    eng = energy_expression(csf; verbosity=verbosity-3)[2][1]
    hf_equations(csf, eng; verbosity=verbosity)
end
=#

function hf_equations(config::Configuration{O}; verbosity=0,
                      selector::Function = peel) where {O<:SpinOrbital}
    verbosity > 0 &&
        println("Deriving HF equations for $(config)")
    h = one_body_hamiltonian_matrix(O, [config], selector=selector)[1]
    HC = two_body_hamiltonian_matrix(O, [config], selector=selector)[1]

    orbitals = selector(config).orbitals
    symmetries = find_symmetries(orbitals)

    map(orbitals) do orb
        corb = Conjugate(orb)

        ∂h = filter(d -> !iszero(d), diff.(h.integrals, Ref(corb)))
        ∂HC = diff.(HC.integrals, Ref(corb))
        nz = .!map(iszero, ∂HC)

        multipole_terms = multipole_expand.(HC.integrals[nz])

        # Find all other orbitals of the same symmetry as the current
        # one. These will be used to create a projector, that projects
        # out their components.
        symmetry_orbitals = filter(o -> o != orb, symmetries[symmetry(orb)])

        orb,(∂h, (∂HC[nz], multipole_terms)),symmetry_orbitals
    end
end

"""
    diff(atom)

Differentiate the energy expression associated with the `atom`'s
CSF(s) with respect to the atomic orbitals to derive the Hartree–Fock
equations for the orbitals.
"""
function Base.diff(atom::Atom; verbosity=0)
    length(atom.configurations) > 1 &&
        throw(ArgumentError("Cannot derive Hartree–Fock equations for a multi-configurational atom"))

    map(hf_equations(atom.configurations[1]; verbosity=verbosity)) do (orb,equation,symmetry_orbitals)
        verbosity > 3 && display(equation)
        hfeq = HFEquation(atom, equation, orb, symmetry_orbitals)
        verbosity > 2 && println(hfeq)
        hfeq
    end
end

Base.view(fock::Fock{A,E}, args...) where {A<:Atom,E} =
    view(fock.quantum_system.radial_orbitals.mul.factors[2], args...)
