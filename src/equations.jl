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
    energy_matrix!(H, hfeqs::AtomicEquations[, which=:energy])

Compute the energy matrix by computing the energy observable and
storing it in `H`. Requires that `hfeqs` has the `:energy` and
`:kinetic_energy` [`Observable`](@ref)s registered (this is the
default).
"""
function SCF.energy_matrix!(H::HM, hfeqs::AtomicEquations,
                            which::Symbol=:total_energy) where {HM<:AbstractMatrix}
    observable = hfeqs.observables[which]
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

# * Setup orbital equations

function get_operator(::FieldFreeOneBodyHamiltonian, atom::Atom,
                      orbital::O, source_orbital::O) where O
    orbital == source_orbital ||
        throw(ArgumentError("FieldFreeOneBodyHamiltonian not pertaining to $(orbital)"))
    AtomicOneBodyHamiltonian(atom, orbital)
end

for (i,HT) in enumerate([:KineticEnergyHamiltonian, :PotentialEnergyHamiltonian])
    @eval begin
        get_operator(::$HT, atom::Atom, orbital::O, ::O) where O =
            AtomicOneBodyHamiltonian(one_body_hamiltonian(Tuple, atom, orbital)[$i],
                                     orbital)
    end
end

function get_operator(op::CoulombPotentialMultipole, atom::Atom,
                      orbital::O, ::O) where O
    a,b = op.a[1],op.b[1]
    HFPotential(:exchange, op.o.k, a, a, view(atom, a), view(atom, a))
end

get_operator(top::IdentityOperator, atom::Atom,
             ::O, source_orbital::O) where O =
                 SourceTerm(top, source_orbital, view(atom, source_orbital))

get_operator(op::QO, atom::Atom, ::O, ::O) where {QO,O} =
    throw(ArgumentError("Unsupported operator $op"))

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

        for (integral,equation_terms) in equation.terms
            if integral > 0
                operator = get_integral(integrals, integral_map, eqs.integrals[integral])
                pushterms!(terms, operator, equation_terms,
                           integrals, integral_map, eqs.integrals)
            else
                for t in equation_terms
                    operator = get_operator(t.operator, atom, orbital, t.source_orbital)

                    pushterms!(terms, operator, [t],
                               integrals, integral_map, eqs.integrals)
                end
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
                   selector::Function=cfg -> outsidecoremodel(cfg, atom.potential),
                   observables::Dict{Symbol,Tuple{<:QuantumOperator,Bool}} = Dict(
                       :total_energy => (H,false),
                       :double_counted_energy => (H,true),
                       :kinetic_energy => (KineticEnergyHamiltonian(),false),
                   ),
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

        if verbosity > 1
            println("Common integrals:")
            display(eqs.integrals)
            println()
        end
    end

    integrals = Vector{OrbitalIntegral# {<:Any,O,T,B,RadialOrbital{T,B}}
                       }()
    integral_map = Dict{Any,Int}()
    append_common_integrals!(integrals, integral_map,
                             atom, eqs.integrals)

    hfeqs = generate_atomic_orbital_equations(atom, eqs,
                                              integrals, integral_map,
                                              symmetries; verbosity=verbosity)

    observables = map(collect(pairs(observables))) do (k,(operator,double_counted))
        k => Observable(operator, atom, overlaps,
                       integrals, integral_map,
                       symmetries, selector,
                       double_counted=double_counted)
    end |> Dict

    AtomicEquations(atom, hfeqs, integrals, observables)
end
