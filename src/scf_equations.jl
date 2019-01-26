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

include("symbolic_hfequations.jl")

function HFEquation(atom::A, (one_body,(two_body,multipole_terms))::E,
                    orbital::O) where {T,ΦT, #<:RadialCoeff{T},
                                       B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       A<:Atom{T,ΦT,B,O},
                                       E<:Tuple{Vector{<:OneBodyHamiltonian},Tuple{Vector{<:DirectExchangePotentials},Vector{NTuple{2,Vector{Pair{Int,Float64}}}}}}}
    R = radial_basis(atom)

    ĥ = one_body_hamiltonian(atom, orbital)

    korb = Ket(orbital)

    OV = typeof(view(atom, 1))
    PO = typeof(R*Diagonal(Vector{T}(undef, size(R,2)))*R')
    direct_potentials = Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}()
    exchange_potentials = Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}()

    count(iszero.(one_body)) == 1 ||
        throw(ArgumentError("There can only be one one-body Hamiltonian per orbital"))

    all(AngularMomentumAlgebra.isdiagonal.(two_body)) ||
        throw(ArgumentError("Non-diagonal repulsion potentials not yet supported"))

    for (tb,mpt) ∈ zip(two_body,multipole_terms)
        tb.o.v == orbital ||
            throw(ArgumentError("Repulsion potential $(tb) not pertaining to $(orbital)"))

        other = tb.a
        v = view(atom, other)

        for (k,c) in mpt[1]
            push!(direct_potentials,
                  HFPotential(:direct, k, other, v) => c)
        end
        for (k,c) in mpt[2]
            push!(exchange_potentials,
                  HFPotential(:exchange, k, other, v) => c)
        end
    end

    hamiltonian = OrbitalSplitHamiltonian(ĥ, direct_potentials, exchange_potentials)

    HFEquation(atom, (one_body,two_body), orbital, view(atom, orbital), hamiltonian)
end

SCF.energy(hfeq::HFEquation{E,O,M}) where {E,O,M} = (hfeq.ϕ' * hfeq.hamiltonian * hfeq.ϕ)[1]

function Base.show(io::IO, hfeq::HFEquation)
    write(io, "Hartree–Fock equation: 0 = [𝓗  - E($(hfeq.orbital))]|$(hfeq.orbital)⟩ = [$(hfeq.hamiltonian) - E($(hfeq.orbital))]|$(hfeq.orbital)⟩")
    EHa = SCF.energy(hfeq)
    write(io, "\n    ⟨$(hfeq.orbital)| 𝓗 |$(hfeq.orbital)⟩ = $(EHa) Ha = $(27.211EHa) eV")
end

# * Setup Hartree–Fock equations

function hf_equations(csf::NonRelativisticCSF, eng::Number; verbosity=0)
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

function hf_equations(config::Configuration{O}; verbosity=0) where {O<:SpinOrbital}
    verbosity > 0 &&
        println("Deriving HF equations for $(config)")
    h = one_body_hamiltonian_matrix(O, [config])[1]
    HC = two_body_hamiltonian_matrix(O, [config])[1]
    map(config.orbitals) do orb
        corb = Conjugate(orb)
        multipole_terms = multipole_expand.(HC.integrals)
        orb,(diff.(h.integrals, Ref(corb)), (diff.(HC.integrals, Ref(corb)), multipole_terms))
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

    map(hf_equations(atom.configurations[1]; verbosity=verbosity)) do (orb,equation)
        verbosity > 3 && display(equation)
        hfeq = HFEquation(atom, equation, orb)
        verbosity > 2 && println(hfeq)
        hfeq
    end
end

Base.view(fock::Fock{A,E}, args...) where {A<:Atom,E} =
    view(fock.quantum_system.radial_orbitals.mul.factors[2], args...)
