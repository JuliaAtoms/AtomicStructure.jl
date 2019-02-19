"""
    Observable

Represents a physical quantity that can be observed, which is
calculated as the matrix element of an operator between two
configurations. All physical observables are real.
"""
mutable struct Observable{T,B,Equations<:AbstractVector{<:AtomicOrbitalEquation}}
    equations::Equations
    tmp::RadialOrbital{T,B}
end

"""
    Observable(operator, atom, overlaps, integrals)

Construct an observable corresponding the `operator` acting on
`atom`. `overlaps` is a list of non-orthogonal, `integrals` a list of
common integrals, and `integral_map` is a mapping from symbolic
integrals to [`OrbitalIntegral`](@ref)s.
"""
function Observable(operator::QuantumOperator, atom::A,
                    overlaps::Vector{<:OrbitalOverlap},
                    integrals::Vector{OrbitalIntegral},
                    integral_map::Dict{Any,Int},
                    symmetries::Dict) where {A<:Atom}
    M = Matrix(operator, atom.configurations, overlaps)

    symbolic_integrals = []

    m,n = size(M)

    equation_system = map(atom.orbitals) do orbital
        corb = conj(orbital)
        equation = orbital_equation(M, corb, symbolic_integrals)

        # We need to filter out those terms that are dependent on
        # corb, to avoid double-counting when computing the
        # observables.
        for i = 1:m
            for j = 1:n
                iszero(M[i,j]) && continue
                M[i,j] = NBodyMatrixElement(filter(t -> !isdependent(t, corb),
                                                   M[i,j].terms))
            end
        end
        equation
    end

    eqs = generate_atomic_orbital_equations(
        atom, MCEquationSystem(equation_system, symbolic_integrals),
        integrals, integral_map,
        symmetries)

    Observable(eqs, similar(first(eqs).ϕ))
end

"""
    observe!(A::M, o::Observable)

Compute the observable `o` between all configurations and store the
results as matrix elements of `A`.
"""
function observe!(A::M, o::Observable{T}) where {M<:AbstractMatrix{<:Real},T}
    A .= zero(T)
    for eq in o.equations
        for term in eq.hamiltonian.terms
            materialize!(MulAdd(coefficient(term), term.A, eq.ϕ, zero(T), o.tmp))
            A[term.i,term.j] += real((eq.ϕ'o.tmp)[1])
        end
    end
    A
end
