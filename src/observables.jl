"""
    Observable

Represents a physical quantity that can be observed, which is
calculated as the matrix element of an operator between two
configurations. All physical observables are real.
"""
mutable struct Observable{T,B<:Basis,Equations<:AbstractVector{<:AtomicOrbitalEquation}}
    equations::Equations
    tmp::RadialOrbital{T,B}
end

"""
    Observable(operator, atom, overlaps, integrals[; double_counted=false])

Construct an observable corresponding the `operator` acting on `atom`;
if `double_counted`, only return those terms that would be
double-counted, otherwise return the normal observable
equations. `overlaps` is a list of non-orthogonal, `integrals` a list
of common integrals, and `integral_map` is a mapping from symbolic
integrals to [`OrbitalIntegral`](@ref)s.
"""
function Observable(operator::QuantumOperator, atom::A,
                    overlaps::Vector{<:OrbitalOverlap},
                    integrals::Vector{OrbitalIntegral},
                    integral_map::Dict{Any,Int},
                    symmetries::Dict,
                    selector::Function;
                    double_counted::Bool=false) where {T,A<:Atom{T}}
    M = Matrix(operator, selector.(atom.configurations), overlaps)

    m,n = size(M)

    if double_counted
        numbodies(operator) > 2 &&
            throw(ArgumentError("Cannot construct double-counted observables for N-body operators with N > 2."))
        M *= one(T)/2
        # We need to filter out those terms that are due to zero- and
        # one-body interactions, since they would not be
        # double-counted.
        for i = 1:m
            for j = 1:n
                iszero(M[i,j]) && continue
                M[i,j] = NBodyMatrixElement(filter(t -> any(f -> numbodies(f) > 1,
                                                           t.factors),
                                                   M[i,j].terms))
            end
        end
    end

    symbolic_integrals = []

    equation_system = map(atom.orbitals) do orbital
        corb = conj(orbital)
        equation = orbital_equation(M, corb, symbolic_integrals)

        # We need to filter out those terms that are dependent on
        # corb, to avoid double-counting when computing the
        # observables.
        for i = 1:m
            for j = 1:n
                (iszero(M[i,j]) || double_counted) && continue
                M[i,j] = NBodyMatrixElement(filter(t -> !isdependent(t, corb),
                                                   M[i,j].terms))
            end
        end
        equation
    end

    # TODO: Should observables really project out symmetry orbitals?
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
            A[term.i,term.j] += real((materialize(applied(*, eq.ϕ', o.tmp))))
        end
    end
    A
end
