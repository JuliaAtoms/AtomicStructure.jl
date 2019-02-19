"""
    push_common_integral!(integrals, integral_map,
                          integral, atom)

Push the `integral` to `integrals`, constructing the correct
[`OrbitalIntegral`](@ref) pertaining to the `atom`. Record the index
of `integral` in `integrals` in the `integral_map`.
"""
function push_common_integral!(integrals::Vector{OrbitalIntegral},
                               integral_map::Dict{Any,Int},
                               integral, atom::A) where {A<:Atom}
    if integral isa OrbitalOverlap
        a,b = integral.a,integral.b
        push!(integrals, OrbitalOverlapIntegral(a, b, view(atom, a), view(atom, b)))
    elseif integral isa CoulombPotentialMultipole
        a,b = integral.a[1],integral.b[1]
        k = integral.o.k
        # All CoulombPotentialMultipoles that occur as common
        # integrals in an equation are direct potentials acting on
        # either the orbital, or on a source orbital coupled
        # through configuration interaction.
        push!(integrals, HFPotential(:direct, k, a, b, view(atom, a), view(atom, b)))
    else
        throw(ArgumentError("$(integral) of type $(typeof(integral)) not yet supported"))
    end
    integral_map[integral] = length(integrals)
end

function append_common_integrals!(integrals::Vector{OrbitalIntegral},
                                  integral_map::Dict{Any,Int},
                                  atom::A, equation_integrals) where {A<:Atom}
    for integral in equation_integrals
        integral ∈ keys(integral_map) && continue
        push_common_integral!(integrals, integral_map, integral, atom)
    end
end

get_integral(integrals::Vector{OrbitalIntegral},
             integral_map::Dict{Any,Int},
             symbolic_integral) =
    integrals[integral_map[symbolic_integral]]

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
function pushterms!(terms::Vector{<:OrbitalHamiltonianTerm{O,T,B,OV}},
                    operator::QO,
                    equation_terms::Vector,
                    integrals::Vector{OrbitalIntegral},
                    integral_map::Dict{Any,Int},
                    symbolic_integrals) where {O,T,B,OV,QO}
    for eq_term in equation_terms
        push!(terms,
              OrbitalHamiltonianTerm{O,T,B,OV,QO}(
                  eq_term.i, eq_term.j, T(eq_term.coeff),
                  operator,
                  [get_integral(integrals, integral_map, symbolic_integrals[i])
                   for i in eq_term.integrals]))
    end
end
