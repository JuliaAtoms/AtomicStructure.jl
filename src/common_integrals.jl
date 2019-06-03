get_integral(integrals::Vector, integral_map::Dict, symbolic_integral) =
    integrals[integral_map[symbolic_integral]]

function get_or_create_integral!(integrals, integral_map, symbolic_integral, atom)
    if symbolic_integral ∈ keys(integral_map)
        get_integral(integrals, integral_map, symbolic_integral)
    else
        i = create_integral(symbolic_integral, atom, integrals, integral_map)
        isnothing(i) && return
        push!(integrals, i)
        integral_map[symbolic_integral] = length(integrals)
        i
    end
end

function create_integral(symbolic_integral::OrbitalOverlap, atom::Atom, integrals, integral_map)
    a,b = symbolic_integral.a,symbolic_integral.b
    OrbitalOverlapIntegral(a, b, view(atom, a), view(atom, b))
end

function create_integral(symbolic_integral::CoulombPotentialMultipole, atom::Atom, integrals, integral_map)
    a,b = symbolic_integral.a[1],symbolic_integral.b[1]
    k = symbolic_integral.o.k
    # All CoulombPotentialMultipoles that occur as common
    # integrals in an equation are direct potentials acting on
    # either the orbital, or on a source orbital coupled
    # through configuration interaction.
    HFPotential(:direct, k, a, b, view(atom, a), view(atom, b))
end

function create_integral(symbolic_integral::OrbitalMatrixElement{2,<:SpinOrbital,CoulombInteractionMultipole,<:SpinOrbital},
                         atom::Atom, integrals, integral_map)
    a,b = symbolic_integral.a[1],symbolic_integral.b[1]
    ∂I = diff(symbolic_integral, conj(a))
    # The operator in a two-electron integral is a HFPotential (that
    # may be shared with many other two-electron integrals). The
    # potential must be updated before the two-electron integral is
    # computed, and we achieve this by making sure it appears before
    # the two-electron integral in the list of integrals (which is
    # updated sequentially).
    hfpotential = get_or_create_integral!(integrals, integral_map, ∂I.operator, atom)
    # We know that the coefficient is purely numeric (and likely
    # unity), since we varied a single OrbitalMatrixElement with
    # respect to one of the orbitals.
    coeff = ∂I.factor.coeff
    OperatorMatrixElement(a, b, view(atom, a), view(atom, b), hfpotential.V̂, coeff)
end

create_integral(symbolic_integral, atom::Atom, integrals, integral_map) =
        throw(ArgumentError("$(symbolic_integral) of type $(typeof(symbolic_integral)) not yet supported"))

function append_common_integrals!(integrals::Vector, integral_map::Dict, atom::Atom, equation_integrals)
    for symbolic_integral in equation_integrals
        get_or_create_integral!(integrals, integral_map, symbolic_integral, atom)
    end
end

"""
    pushterms!(terms, operator, equation_terms,
               integrals, integral_map, symbolic_integrals)

For each term in `equation_terms`, push a term, located at CI
coordinates `i,j`, of the overall orbital Hamiltonian to `terms`,
constructed from `operator` and a product of orbital `integrals`,
multiplied by an overall factor given by expression and multipole
expansions. `integrals` contain common [`OrbitalIntegral`](@ref)s
and `integral_map` maps from `symbolic_integrals` to `integrals`.
"""
function pushterms!(terms::Vector{<:OrbitalHamiltonianTerm{O,O,T,B,OV}},
                    operator::QO,
                    equation_terms::Vector,
                    integrals::Vector{OrbitalIntegral},
                    integral_map::Dict{Any,Int},
                    symbolic_integrals) where {O,T,B,OV,QO}
    for eq_term in equation_terms
        push!(terms,
              OrbitalHamiltonianTerm{O,O,T,B,OV,QO}(
                  eq_term.i, eq_term.j, T(eq_term.coeff),
                  operator,
                  [get_integral(integrals, integral_map, symbolic_integrals[i])
                   for i in eq_term.integrals]))
    end
end
