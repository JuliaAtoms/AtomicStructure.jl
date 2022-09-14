function evaluate(oo::OrbitalOverlap, atom::Atom{T}, integrals, integral_map) where {N,T}
    a,b = oo.a,oo.b
    av,bv = view(atom,a),view(atom,b)
    dot(av.args[2], atom.S̃, bv.args[2])
end

function evaluate(ome::OrbitalMatrixElement{N}, atom::Atom{T}, integrals, integral_map) where {N,T}
    a,o,b = ome.a,ome.o,ome.b
    op = if N == 1
        get_operator(o, atom, a[1], b[1])
    elseif N == 2
        symbolic_integral = diff(ome, conj(a[1])).operator
        get_or_create_integral!(integrals, integral_map, symbolic_integral, atom)
    end
    iszero(op) && return zero(T)
    av,bv = view(atom,a[1]),view(atom,b[1])
    tmp = similar(bv)
    tmp.args[2] .= zero(T)
    materialize!(MulAdd(one(T), op, bv, zero(T), tmp))
    dot(av.args[2], atom.S̃, tmp.args[2])
end

evaluate(nbt::NBodyTerm, atom::Atom, integrals, integral_map) =
    nbt.coeff*prod(evaluate(f, atom, integrals, integral_map) for f in nbt.factors)

evaluate(nbme::NBodyMatrixElement, atom::Atom, integrals, integral_map) =
    sum(evaluate(t, atom, integrals, integral_map) for t in nbme.terms)

function evaluate_matrix!(H::AbstractMatrix, E::EnergyExpression, atom::Atom)
    size(H) == size(E) || throw(DimensionMismatch("Matrix sizes do not agree"))

    eqs = diff(E, conj.(atom.orbitals))
    integrals, integral_map, poisson_cache = common_integrals(atom, eqs.integrals)

    rows = rowvals(E)
    vals = nonzeros(E)
    for j = 1:size(E,2)
        for k in nzrange(E, j)
            i = rows[k]
            H[i,j] = evaluate(vals[k], atom, integrals, integral_map)
        end
    end
    H
end

function evaluate_matrix(E::EnergyExpression, atom::Atom)
    m,n = size(E)
    evaluate_matrix!(zeros(eltype(atom), m, n), E, atom)
end

evaluate_matrix(O::QuantumOperator, atom::Atom; kwargs...) =
    evaluate_matrix(Matrix(O, atom; kwargs...), atom)

energy_matrix(atom::Atom; kwargs...) =
    evaluate_matrix(atomic_hamiltonian(atom), atom; kwargs...)

function expectation_value(O, atom::Atom; kwargs...)
    M = evaluate_matrix(O, atom; kwargs...)
    c = SCF.coefficients(atom)
    dot(c, M, c)
end

energy(atom::Atom; kwargs...) =
    expectation_value(atomic_hamiltonian(atom), atom; kwargs...)

export evaluate_matrix!, evaluate_matrix, energy_matrix,
    expectation_value, energy
