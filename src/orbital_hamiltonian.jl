# * Split Hamiltonian

"""
    OrbitalHamiltonianTerm(i, j, coeff, A, integrals)

Represents a term in the orbital Hamiltonian arising from a variation
of the energy expressions between configurations `i` and `j` in the
multi-configurational expansion. `coeff` is the numeric coefficient,
`A` is the operator acting on the orbital, and `integrals` is a vector
of [`OrbitalIntegral`](@ref)s arising from the presence of
non-orthogonal orbitals and whose values should be multiplied to form
the overall coefficient.
"""
struct OrbitalHamiltonianTerm{aO,bO,T,B<:Basis,OV,QO}
    i::Int
    j::Int
    coeff::T
    A::QO
    integrals::Vector{OrbitalIntegral{<:Any,aO,bO,T,B,OV}}
end
"""
    coefficient(term::OrbitalHamiltonianTerm)

Return the multiplicative coefficient pertaining to `term`,
_excluding_ the `conj(c_i)*c_j` mixing coefficients, due to the
configuration-interaction.
"""
coefficient(term::OrbitalHamiltonianTerm) =
    term.coeff*prod(integral_value.(term.integrals))

"""
    coefficient(term::OrbitalHamiltonianTerm, c::Vector)

Return the multiplicative coefficient pertaining to `term`,
_including_ the `conj(c_i)*c_j` mixing coefficients, due to the
configuration-interaction.
"""
function coefficient(term::OrbitalHamiltonianTerm, c::Vector)
    coeff = coefficient(term)
    # We use zero-valued i,j coordinates to designate terms which do
    # not originate from the configuration expansion, but rather added
    # in manually, such as energy shifts.
    (iszero(term.i) || iszero(term.j)) && return coeff

    conj(c[term.i])*c[term.j]*coeff
end

"""
    OrbitalHamiltonian(R, terms, mix_coeffs, projector, orbital)

The Hamiltonian for `orbital` is constructed from a radial basis `R`,
a set of [`OrbitalHamiltonianTerm`](@ref) `terms` that describe the
various interactions between orbitals, `mix_coeffs` which are the
mixing coefficents for the multi-configurational expansion. The
`projector` ensures orthogonality between orbital pairs which have
Lagrange multipliers associated with them, by projecting out
components of other orbitals every time the `OrbitalHamiltonian`
action on `orbital` is computed.
"""
mutable struct OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj,RT}
    R::RT
    terms::Vector{OrbitalHamiltonianTerm{aO,bO,T,B,OV}}
    mix_coeffs::Vector{T}
    projector::Proj
    orbital::O
end

Base.axes(hamiltonian::OrbitalHamiltonian) =
    (axes(hamiltonian.R,1),axes(hamiltonian.R,1))

Base.axes(hamiltonian::OrbitalHamiltonian, i) =
    axes(hamiltonian)[i]

Base.eltype(hamiltonian::OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj}) where {aO,bO,O,T,B,OV,Proj} = T

"""
    energy_matrix!(H, hamiltonian, ϕ)

Compute the contribution of `hamiltonian` to the Hamiltonian matrix
`H` by repeatedly acting on the associated radial orbital `ϕ` with the
different multi-configurational [`OrbitalHamiltonianTerm`](@ref)s of
`hamiltonian`.
"""
function SCF.energy_matrix!(H::HM, hamiltonian::OrbitalHamiltonian{aO,bO,O,T,B},
                            ϕ::RadialOrbital{T,B}) where {HM<:AbstractMatrix,aO,bO,O,T,B}
    tmp = similar(ϕ)
    for term in hamiltonian.terms
        materialize!(MulAdd(coefficient(term), term.A, ϕ, zero(T), tmp))
        H[term.i,term.j] += (ϕ'tmp)[1]
    end
    H
end

"""
    filter(fun::Function, H::OrbitalHamiltonian)

Filter the [`OrbitalHamiltonianTerm`](@ref)s of `H` according to the
predicate `fun`.
"""
Base.filter(fun::Function, H::OrbitalHamiltonian) =
    OrbitalHamiltonian(H.R, filter(fun, H.terms),
                       H.mix_coeffs, H.projector, H.orbital)

function Base.getindex(H::OrbitalHamiltonian, term::Symbol)
    term == :total && return H
    operator_type = if term == :onebody
        AtomicOneBodyHamiltonian
    elseif term == :direct
        DirectPotential
    elseif term == :exchange
        ExchangePotential
    elseif term == :source
        SourceTerm
    elseif term == :shift
        ShiftTerm
    else
        throw(ArgumentError("Unknown Hamiltonian term $term"))
    end

    filter(t -> t.A isa operator_type, H)
end

# ** Materialization

const OrbitalHamiltonianMatrixElement{aO,bO,O,T,B<:Basis} =
    Mul{<:Any,<:Tuple{<:AdjointRadialOrbital{T,B},
                      <:OrbitalHamiltonian{aO,bO,O,T,B},
                      <:RadialOrbital{T,B}}}

const OrbitalHamiltonianMatrixVectorProduct{aO,bO,O,T,B<:Basis} =
    Mul{<:Any,<:Tuple{<:OrbitalHamiltonian{aO,bO,O,T,B},<:RadialOrbital{T,B}}}

const OrbitalHamiltonianMatrixMatrixProduct{aO,bO,O,T,B<:Basis} =
    Mul{<:Any,<:Tuple{<:OrbitalHamiltonian{aO,bO,O,T,B},<:RadialOrbitals{T,B}}}

Base.eltype(::OrbitalHamiltonianMatrixVectorProduct{aO,bO,O,T}) where {aO,bO,O,T} = T

function Base.copyto!(dest::RadialOrbital{T,B},
                      matvec::OrbitalHamiltonianMatrixVectorProduct{aO,bO,O,T,B}) where {aO,bO,O,T,B<:Basis}
    axes(dest) == axes(matvec) || throw(DimensionMismatch("axes must be the same"))
    R′,v = dest.args
    hamiltonian,b = matvec.args

    c = hamiltonian.mix_coeffs

    v .= zero(T)

    # Compute the additive action of each of the terms in the
    # Hamiltonian, weighting them by their coefficients which are
    # constructed from
    # 1) the configuration-interaction mixing coefficients,
    # 2) coefficients due to the energy expression, and
    # 3) values of various orbitals integrals resulting from
    #    non-orthogonal orbitals.
    for term in hamiltonian.terms
        coeff = coefficient(term, c)
        materialize!(MulAdd(coeff, term.A, b, one(T), dest))
    end

    # Project out all components parallel to other orbitals of the
    # same symmetry.
    projectout!(dest, hamiltonian.projector)

    dest
end

function Base.copyto!(dest::RadialOrbitals{T,B},
                      matvec::OrbitalHamiltonianMatrixMatrixProduct{aO,bO,O,T,B}) where {aO,bO,O,T,B<:Basis}
    axes(dest) == axes(matvec) || throw(DimensionMismatch("axes must be the same"))
    R′,dv = dest.args
    hamiltonian,b = matvec.args
    bv = b.args[2]

    n = size(dv,2)
    for j = 1:n
        copyto!(applied(*, hamiltonian.R, view(dv, :, j)),
                hamiltonian ⋆ applied(*, hamiltonian.R, view(bv, :, j)))
    end
    dest
end

Base.similar(matvec::OrbitalHamiltonianMatrixVectorProduct{aO,bO,O,T,B}) where {aO,bO,O,T,B<:Basis} =
    similar(matvec.args[2])

LazyArrays.materialize(matvec::OrbitalHamiltonianMatrixVectorProduct{aO,bO,O,T,B}) where {aO,bO,O,T,B<:Basis,V} =
    copyto!(similar(matvec), matvec)

function LazyArrays.materialize(matel::OrbitalHamiltonianMatrixElement{aO,bO,O,T,B}) where {aO,bO,O,T,B<:Basis}
    a,op,b = matel.args
    materialize(applied(*, a, materialize(op⋆b)))
end

# *** Materialization into a matrix

"""
    copyto!(dest::AbstractMatix, hamiltonian::OrbitalHamiltonian)

Materialize the orbital `hamiltonian` into matrix form and store it in
`dest`, using the current values of all other orbitals. This is only
possible if the orbital `hamiltonian` does *not* contain any
[`ExchangePotential`](@ref)s or [`SourceTerm`](@ref)s, since the
former is non-local (and thus not representable as a matrix) and the
latter is not a linear operator (but an affine one).

Typical usage is to compute an easily factorizable matrix that can be
used for preconditioning the solution of the full equation.
"""
function Base.copyto!(dest::M, hamiltonian::OrbitalHamiltonian) where {T,M<:AbstractMatrix{T}}
    m = size(hamiltonian.R,2)
    size(dest) == (m,m) || throw(DimensionMismatch("axes must be the same"))

    c = hamiltonian.mix_coeffs

    dest .= zero(T)

    for term in hamiltonian.terms
        (term.A isa SourceTerm ||
         term.A isa ExchangePotential) &&
         throw(ArgumentError("It is not possible to materialize a $(typeof(term.A)) as a matrix"))
        coeff = coefficient(term, c)

        op = if term.A isa AtomicOneBodyHamiltonian
            term.A.op.args[2]
        elseif term.A isa DirectPotential
            term.A.V̂.args[2]
        elseif term.A isa ShiftTerm
            term.A.shift
        else
            throw(ArgumentError("Don't know how to materialize a $(typeof(term.A)) as a matrix"))
        end

        dest += coeff*op
    end

    dest
end

function Base.similar(h::OrbitalHamiltonian{aO,bO,O,T,B}, ::Type{T}) where {aO,bO,O,T,B<:AbstractFiniteDifferences}
    R = h.R
    m = size(R,2)
    # TODO: This is only valid for RadialDifferences of
    # FiniteDifferencesQuasi.jl
    o = ones(T, m)
    SymTridiagonal(o,0*o[2:end])
end

Base.similar(h::OrbitalHamiltonian{aO,bO,O,T,B}, ::Type{T}) where {aO,bO,O,T,B<:BasisOrRestricted{<:FEDVR}} =
    Matrix(undef, h.R)

LazyArrays.materialize(h::OrbitalHamiltonian) =
    copyto!(similar(h, eltype(h)), h)

# ** Arithmetic

"""
    h::OrbitalHamiltonian + λ::UniformScaling

Shift the [`OrbitalHamiltonian`](@ref) `h` by `λ`.
"""
function Base.:(+)(h::OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj,RT}, λ::UniformScaling) where {aO,bO,O,T,B,OV,Proj,RT}
    # The zeros designate that the shift is not to be weigthed by the mixing coefficients
    shift_term = OrbitalHamiltonianTerm(0, 0, one(T), ShiftTerm(λ),
                                        Vector{OrbitalIntegral{<:Any,aO,bO,T,B,OV}}())
    OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj,RT}(h.R, vcat(h.terms, shift_term),
                                         h.mix_coeffs, h.projector, h.orbital)
end

"""
    h::OrbitalHamiltonian - λ::UniformScaling

Shift the [`OrbitalHamiltonian`](@ref) `h` by `-λ`.
"""
Base.:(-)(h::OrbitalHamiltonian, λ::UniformScaling) = h + (-λ)

# ** Krylov wrapper

"""
    SCF.KrylovWrapper(hamiltonian::OrbitalHamiltonian)

Construct a `KrylovWrapper` such that `hamiltonian`, that acts on
function spaces, can be used in a Krylov solver, which works with
linear algebra vector spaces.
"""
SCF.KrylovWrapper(hamiltonian::OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj}) where {aO,bO,O,T,B,OV,Proj} =
    KrylovWrapper{T,OrbitalHamiltonian{aO,bO,O,T,B,OV,Proj}}(hamiltonian)

Base.size(hamiltonian::OrbitalHamiltonian, ::SCF.KrylovWrapper) =
    (size(hamiltonian.R,2),size(hamiltonian.R,2))

"""
    mul!(y, A::KrylovWrapper{T,<:OrbitalHamiltonian}, x)

Materialize the action of the [`OrbitalHamiltonian`](@ref) on the
linear algebra vector `x` and store the result in `y`, by wrapping
them both with the `QuasiMatrix` necessary to transform `x` and `y` to
the function space of the Hamiltonian.
"""
LinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,Hamiltonian}, x::V₂) where {V₁,V₂,T,B,Hamiltonian<:OrbitalHamiltonian} =
    copyto!(applied(*, A.hamiltonian.R, y),
            A.hamiltonian⋆(applied(*, A.hamiltonian.R, x)))

# ** Preconditioner

"""
    IterativeFactorizations.preconditioner(hamiltonian::OrbitalHamiltonian)

Return a factorization of the matrix corresponding to `hamiltonian`,
where all terms arising from exchange and configuration interaction
have been removes, since they cannot be represented by a matrix.
"""
function IterativeFactorizations.preconditioner(hamiltonian::OrbitalHamiltonian)
    # To form the preconditioner, we select all terms of the shifted
    # Hamiltonian, except the exchange potentials and source terms,
    # which are not factorizable.
    Ph = filter(t -> !(t.A isa ExchangePotential || t.A isa SourceTerm), hamiltonian)
    hm = materialize(Ph)
    factorize(hm)
end
