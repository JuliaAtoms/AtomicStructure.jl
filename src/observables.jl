"""
    Observable

Represents a physical quantity that can be observed, which is
calculated as the matrix element of an operator between two
configurations. All physical observables are real.
"""
mutable struct Observable{T,B<:Basis,Equations<:AbstractVector{<:AtomicOrbitalEquation},Metric}
    equations::Equations
    tmp::RadialOrbital{T,B}
    S̃::Metric
end

"""
    Observable(operator, atom, overlaps, integrals, integral_map, selector
               [; double_counted=false])

Construct an observable corresponding the `operator` acting on `atom`;
if `double_counted`, only return those terms that would be
double-counted, otherwise return the normal observable
equations. `overlaps` is a list of non-orthogonal, `integrals` a list
of common integrals, and `integral_map` is a mapping from symbolic
integrals to [`OrbitalIntegral`](@ref)s. `selector` selects those
orbitals not modelled by the potential, e.g. all orbitals in case of a
nuclear potential, fewer in case of a
pseudopotential. `double_counted` can be used to only return those
terms of the sum that would be double-counted if simply summing over
orbital contributions (applicable to two-body operators only). Half
the result of the double-counted term is then subtracted the sum over
orbital contributions; this makes the expression slightly more
symmetric in orbital space than if the sum is derived avoiding
double-counting.
"""
function Observable(operator::QuantumOperator, atom::A,
                    overlaps::Vector{<:OrbitalOverlap},
                    integrals::Vector,
                    integral_map::Dict{Any,Int},
                    symmetries::Dict,
                    selector::Function;
                    double_counted::Bool=false,
                    kwargs...) where {T,A<:Atom{T}}
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
        symmetries; kwargs...)

    tmp = similar(first(eqs).ϕ)
    tmp.args[2] .= false # To clear any NaNs

    Observable(eqs, tmp, atom.S̃)
end

realifreal(::Type{R}, v) where {R<:Real} = real(v)
realifreal(::Type{R}, v) where {R} = v

"""
    observe!(A::M, o::Observable)

Compute the observable `o` between all configurations and store the
results as matrix elements of `A`.
"""
function observe!(A::M, o::Observable{T}) where {R,M<:AbstractMatrix{R},T}
    A .= zero(T)
    for eq in o.equations
        for term in eq.hamiltonian.terms
            materialize!(MulAdd(coefficient(term), term.A, eq.ϕ, zero(T), o.tmp))
            A[term.i,term.j] += realifreal(R, dot(eq.ϕ.args[2], o.S̃, o.tmp.args[2]))
        end
    end
    A
end

"""
    observe!(A::M, o::Observable, atom::Atom)

Compute the observable `o` between all configurations, store the
results as matrix elements of `A`, and finally contract with respect
to the mixing coefficients of `atom` and return the result as a
scalar.
"""
function observe!(A::M, o::Observable{T}, atom::Atom) where {R,M<:AbstractMatrix{R},T}
    observe!(A, o)
    c = atom.mix_coeffs
    dot(c, A, c)
end
