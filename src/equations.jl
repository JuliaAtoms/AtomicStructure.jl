# * Atomic System of Equations

"""
    AtomicEquations(atom, equations, integrals)

Structure representing the (e.g. Hartree–Fock) `equations` for `atom`,
along with all `integrals` that are shared between the `equations`.
"""
mutable struct AtomicEquations{T, B<:AbstractQuasiMatrix,
                               O<:AbstractOrbital,
                               A<:Atom{T,B,O},
                               Equations<:AbstractVector{<:AtomicOrbitalEquation}}
    atom::A
    equations::Equations
    integrals::Vector{OrbitalIntegral}
    observables::Dict{Symbol,<:Observable}
end

# Iteration interface
Base.length(hfeqs::AtomicEquations) = length(hfeqs.equations)
Base.iterate(iter::AtomicEquations, args...) = iterate(iter.equations, args...)

"""
    update!(equations::AtomicEquations)

Recompute all integrals using the current values for the radial orbitals.
"""
SCF.update!(equations::AtomicEquations; kwargs...) =
    foreach(integral -> SCF.update!(integral; kwargs...),
            equations.integrals)

"""
    energy_matrix!(H, hfeqs::AtomicEquations)

Compute the energy matrix by computing the energy observable and
storing it in `H`. Requires that `hfeqs` has the `:energy`
[`Observable`](@ref) registered (default).
"""
function SCF.energy_matrix!(H::HM, hfeqs::AtomicEquations) where {HM<:AbstractMatrix}
    observable = hfeqs.observables[:energy]
    observe!(H, observable)
    H
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

"""
    generate_atomic_orbital_equations(atom::Atom, eqs::MCEquationSystem,
                                      integrals, integral_map)

For each variationally derived orbital equation in `eqs`, generate the
corresponding [`AtomicOrbitalEquation`](@ref).
"""
function generate_atomic_orbital_equations(atom::Atom{T,B,O}, eqs::MCEquationSystem,
                                           integrals::Vector{OrbitalIntegral},
                                           integral_map::Dict{Any,Int},
                                           symmetries::Dict;
                                           verbosity=0) where {T,B,O}
    map(eqs.equations) do equation
        orbital = equation.orbital
        terms = Vector{OrbitalHamiltonianTerm{O,T,B,RadialOrbital{T,B}}}()

        ĥ₀ = AtomicOneBodyHamiltonian(atom, orbital)
        pushterms!(terms, ĥ₀, equation.one_body, integrals, integral_map, eqs.integrals)

        for dt in equation.direct_terms
            pushterms!(terms,
                       get_integral(integrals, integral_map, eqs.integrals[dt[1]]),
                       dt[2],
                       integrals, integral_map, eqs.integrals)
        end

        for xt in equation.exchange_terms
            xop = xt[1] # CoulombPotentialMultipole
            a,b = xop.a[1],xop.b[1]
            xoperator = HFPotential(:exchange, xop.o.k, a, a,
                                    view(atom, a), view(atom, a))
            pushterms!(terms, xoperator, [xt[2]], integrals, integral_map, eqs.integrals)
        end

        for st in equation.source_terms
            sop = integrals[st[1]]
            for (coeff,source_orbital) in st[2]
                pushterms!(terms,
                           SourceTerm(sop, source_orbital, view(atom, source_orbital)),
                           [coeff], integrals, integral_map, eqs.integrals)
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

        AtomicOrbitalEquation(atom, equation, orbital, terms, symmetry_orbitals)
    end
end

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
                   H::QuantumOperator=FieldFreeOneBodyHamiltonian()+CoulombInteraction();
                   overlaps::Vector{OrbitalOverlap}=OrbitalOverlap[],
                   selector=cfg -> outsidecoremodel(cfg, atom.potential),
                   observables::Dict{Symbol,<:QuantumOperator} = Dict(:energy => H),
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
    integral_map = Dict{Any,Int}()
    append_common_integrals!(integrals, integral_map,
                             atom, eqs.integrals)

    hfeqs = generate_atomic_orbital_equations(atom, eqs,
                                              integrals, integral_map,
                                              symmetries; verbosity=verbosity)

    observables = map(collect(pairs(observables))) do (k,operator)
        k => Observable(operator, atom, overlaps,
                       integrals, integral_map,
                       symmetries)
    end |> Dict

    AtomicEquations(atom, hfeqs, integrals, observables)
end

export energy
