# * Hartree–Fock equations
# ** Hartree-Fock orbital equation

mutable struct HFEquation{T, B<:AbstractQuasiMatrix,
                          O<:AbstractOrbital,
                          A<:Atom{T,B,O},
                          Equation,
                          M# <:OrbitalHamiltonian{T,B,O}
                          }
    atom::A
    equation::Equation
    orbital::O
    ϕ::RadialOrbital{T,B}
    hamiltonian::M
end

SCF.hamiltonian(hfeq::HFEquation) = hfeq.hamiltonian

function HFEquation(atom::A, equation::Equation, orbital::O,
                    terms::Vector{<:OrbitalHamiltonianTerm{O,T,B}},
                    symmetry_orbitals::Vector{O}) where {T,
                                                         B<:AbstractQuasiMatrix,
                                                         O<:AbstractOrbital,
                                                         A<:Atom{T,B,O},
                                                         Equation}
    OV = typeof(view(atom,orbital))
    projector = Projector(OV[view(atom, other)
                             for other in symmetry_orbitals])

    hamiltonian = OrbitalHamiltonian(radial_basis(atom), terms,
                                     atom.mix_coeffs, projector, orbital)

    HFEquation(atom, equation, orbital, view(atom, orbital), hamiltonian)
end

SCF.energy(hfeq::HFEquation, term=:all) =
    (hfeq.ϕ' * hfeq.hamiltonian[term] * hfeq.ϕ)[1]

SCF.energy_matrix!(H::HM, hfeq::HFEquation) where {HM<:AbstractMatrix} =
    SCF.energy_matrix!(H, hfeq.hamiltonian, hfeq.ϕ)

function Base.show(io::IO, hfeq::HFEquation)
    EHa = SCF.energy(hfeq)
    write(io, "⟨$(hfeq.orbital)| 𝓗 |$(hfeq.orbital)⟩ = $(EHa) Ha = $(27.211EHa) eV")
end

function Base.show(io::IO, ::MIME"text/plain", hfeq::HFEquation)
    write(io, "Hartree–Fock equation: E|$(hfeq.orbital)⟩ = ")
    show(io, hfeq.equation)
    write(io, "\n")
    show(io, hfeq)
    write(io, "\n")
end

# ** Hartree–Fock system of equations

"""
    HFEquations(atom, equations, integrals)

Structure representing the Hartree–Fock `equations` for `atom`, along
with all `integrals` that are shared between the `equations`.
"""
mutable struct HFEquations{T, B<:AbstractQuasiMatrix,
                           O<:AbstractOrbital,
                           A<:Atom{T,B,O},
                           Equations<:AbstractVector{<:HFEquation}}
    atom::A
    equations::Equations
    integrals::Vector{OrbitalIntegral}
end

# Iteration interface
Base.length(hfeqs::HFEquations) = length(hfeqs.equations)
Base.iterate(iter::HFEquations, args...) = iterate(iter.equations, args...)

"""
    update!(equations::HFEquations)

Recompute all integrals using the current values for the radial orbitals.
"""
SCF.update!(equations::HFEquations; kwargs...) =
    foreach(integral -> SCF.update!(integral; kwargs...),
            equations.integrals)

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
"""
    diff(atom[, H]; overlaps=[], selector=outsidecoremodel, verbosity=0)

Differentiate the energy expression of the Hamiltonian `H` associated
with the `atom`'s configurations(s) with respect to the atomic
orbitals to derive the Hartree–Fock equations for the orbitals.

By default, the Hamiltonian
`H=FieldFreeOneBodyHamiltonian()+CoulombInteraction()`.

Non-orthogonality between orbitals can be specified by providing
`OrbitalOverlap`s between these pairs. Only those electrons not
modelled by `atom.potential` of each configuration are considered for
generating the energy expression, this can be changed by choosing
another value for `selector`.
"""
function Base.diff(atom::Atom{T,B,O},
                   H=FieldFreeOneBodyHamiltonian()+CoulombInteraction();
                   overlaps::Vector{OrbitalOverlap}=OrbitalOverlap[],
                   selector=cfg -> outsidecoremodel(cfg, atom.potential),
                   verbosity=0) where {T,B,O}
    configurations = selector.(atom.configurations)
    orbitals = unique_orbitals(configurations)

    energy_expression = Matrix(H, configurations, overlaps)
    symmetries = find_symmetries(orbitals)
    eqs = diff(energy_expression, conj.(orbitals))

    if verbosity > 0
        println("Energy expression:")
        display(energy_expression)
        println()

        println("Hartree–Fock equations:")
        display(eqs)
        println()

        println("Symmetries:")
        display(symmetries)
        println()
    end

    integrals = Vector{OrbitalIntegral# {<:Any,O,T,B,RadialOrbital{T,B}}
                       }()
    for integral in eqs.integrals
        if integral isa OrbitalOverlap
            a,b = integral.a,integral.b
            push!(integrals, OrbitalOverlapIntegral(a, b, view(atom, a), view(atom, b)))
        elseif integral isa CoulombPotentialMultipole
            a,b = integral.a[1],integral.b[1]
            k = integral.o.k
            # All CoulombPotentialMultipoles that occur as common
            # integrals in an equation are direct potentials acting on
            # either the orbital corresponding to the orbital, or on a
            # source orbital coupled through configuration interaction.
            push!(integrals, HFPotential(:direct, k, a, b, view(atom, a), view(atom, b)))
        else
            throw(ArgumentError("$(integral) of type $(typeof(integral)) not yet supported"))
        end
    end

    function pushterms!(terms::Vector{<:OrbitalHamiltonianTerm{O,T,B,OV}},
                        operator::QO,
                        equation_terms::Vector) where {O,T,B,OV,QO}
        for eq_term in equation_terms
            push!(terms,
                  OrbitalHamiltonianTerm{O,T,B,OV,QO}(
                      eq_term.i, eq_term.j, T(eq_term.coeff),
                      operator,
                      [integrals[i] for i in eq_term.integrals]))
        end
    end

    hfeqs = map(eqs.equations) do equation
        orbital = equation.orbital
        ĥ₀ = AtomicOneBodyHamiltonian(atom, orbital)
        terms = Vector{OrbitalHamiltonianTerm{O,T,B,RadialOrbital{T,B}}}()

        pushterms!(terms, ĥ₀, equation.one_body)
        for dt in equation.direct_terms
            pushterms!(terms, integrals[dt[1]], dt[2])
        end
        for xt in equation.exchange_terms
            xop = xt[1] # CoulombPotentialMultipole
            a,b = xop.a[1],xop.b[1]
            xoperator = HFPotential(:exchange, xop.o.k, a, a,
                                    view(atom, a), view(atom, a))
            pushterms!(terms, xoperator, [xt[2]])
        end
        for st in equation.source_terms
            sop = integrals[st[1]]
            for (coeff,source_orbital) in st[2]
                pushterms!(terms, SourceTerm(sop, source_orbital, view(atom, source_orbital)), [coeff])
            end
        end

        # Find all other orbitals of the same symmetry as the current
        # one. These will be used to create a projector, that projects
        # out their components.
        #
        # TODO: Think of what the non-orthogonalities due to
        # `overlaps` imply for the Lagrange multipliers/projectors.
        symmetry_orbitals = filter(!isequal(orbital), symmetries[symmetry(orbital)])
        verbosity > 1 && println("Symmetry: ", symmetry_orbitals)

        HFEquation(atom, equation, orbital, terms, symmetry_orbitals)
    end

    HFEquations(atom, hfeqs, integrals)
end
