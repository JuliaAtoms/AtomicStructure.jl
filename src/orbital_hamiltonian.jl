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
struct OrbitalHamiltonianTerm{O,T,B,OV,QO}
    i::Int
    j::Int
    coeff::T
    A::QO
    integrals::Vector{OrbitalIntegral{<:Any,O,T,B,OV}}
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
mutable struct OrbitalHamiltonian{O,T,B,OV,Proj}
    R::B
    terms::Vector{OrbitalHamiltonianTerm{O,T,B,OV}}
    mix_coeffs::Vector{T}
    projector::Proj
    orbital::O
end

Base.axes(hamiltonian::OrbitalHamiltonian) =
    (axes(hamiltonian.R,1),axes(hamiltonian.R,1))

Base.axes(hamiltonian::OrbitalHamiltonian, i) =
    axes(hamiltonian)[i]

"""
    energy_matrix!(H, hamiltonian, ϕ)

Compute the contribution of `hamiltonian` to the Hamiltonian matrix
`H` by repeatedly acting on the associated radial orbital `ϕ` with the
different multi-configurational [`OrbitalHamiltonianTerm`](@ref)s of
`hamiltonian`.
"""
function SCF.energy_matrix!(H::HM, hamiltonian::OrbitalHamiltonian{O,T,B},
                            ϕ::RadialOrbital{T,B}) where {HM<:AbstractMatrix,O,T,B}
    tmp = similar(ϕ)
    for term in hamiltonian.terms
        materialize!(MulAdd(coefficient(term), term.A, ϕ, zero(T), tmp))
        H[term.i,term.j] += (ϕ'tmp)[1]
    end
    H
end

# function Base.show(io::IO, hamiltonian::OrbitalHamiltonian{T}) where T
#     if iszero(hamiltonian)
#         write(io, "0")
#         return
#     end
#     multiple_terms = sum([!iszero(hamiltonian.ĥ), !isempty(hamiltonian.direct_potentials),
#                           !isempty(hamiltonian.exchange_potentials)]) > 1
#     multiple_terms && write(io, "[")
#     !iszero(hamiltonian.ĥ) && write(io, "ĥ")
#     for (p,c) in hamiltonian.direct_potentials
#         s = sign(c)
#         write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
#     end
#     for (p,c) in hamiltonian.exchange_potentials
#         s = sign(c)
#         write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
#     end
#     multiple_terms && write(io, "]")
#     write(io, "|")
#     show(io, hamiltonian.orbital)
#     write(io, "⟩")
# end

function Base.getindex(H::OrbitalHamiltonian, term::Symbol)
    term == :all && return H
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

    OrbitalHamiltonian(H.R, filter(t -> t.A isa operator_type, H.terms),
                       H.mix_coeffs, H.projector, H.orbital)
end

# ** Materialization

const OrbitalHamiltonianMatrixElement{O,T,B<:AbstractQuasiMatrix,V<:AbstractVector} =
    Mul{<:Tuple,<:Tuple{<:Adjoint{T,V},<:QuasiAdjoint{T,B},<:OrbitalHamiltonian{O,T,B},B,V}}

const OrbitalHamiltonianMatrixVectorProduct{O,T,B<:AbstractQuasiMatrix,V<:AbstractVector} =
    Mul{<:Tuple,<:Tuple{<:OrbitalHamiltonian{O,T,B},B,V}}

Base.eltype(::OrbitalHamiltonianMatrixVectorProduct{O,T}) where {O,T} = T

function Base.copyto!(dest::RadialOrbital{T,B},
                      matvec::OrbitalHamiltonianMatrixVectorProduct{O,T,B,V}) where {O,T,B,V}
    axes(dest) == axes(matvec) || throw(DimensionMismatch("axes must be the same"))
    R′,v = dest.mul.factors
    hamiltonian,R,b = matvec.factors

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
        coeff = conj(c[term.i])*c[term.j]*coefficient(term)
        materialize!(MulAdd(coeff, term.A, R*b, one(T), dest))
    end

    # Project out all components parallel to other orbitals of the
    # same symmetry.
    projectout!(dest, hamiltonian.projector)

    dest
end

function Base.similar(matvec::OrbitalHamiltonianMatrixVectorProduct, ::Type{T}) where T
    A,R,b = matvec.factors
    v = similar(b, T)
    R*v
end

LazyArrays.materialize(matvec::OrbitalHamiltonianMatrixVectorProduct{O,T,B,V}) where {O,T,B,V} =
    copyto!(similar(matvec, eltype(matvec)), matvec)

function LazyArrays.materialize(matel::OrbitalHamiltonianMatrixElement{O,T,B,V}) where {O,T,B,V}
    a,R′,op,R,b = matel.factors
    a*R′*materialize(op⋆R⋆b)
end

# ** Arithmetic

"""
    h::OrbitalHamiltonian + λ::UniformScaling

Shift the [`OrbitalHamiltonian`](@ref) `h` by `λ`.
"""
function Base.:(+)(h::OrbitalHamiltonian{O,T,B,OV,Proj}, λ::UniformScaling) where {O,T,B,OV,Proj}
    # The zeros designate that the shift is not to be weigthed by the mixing coefficients
    shift_term = OrbitalHamiltonianTerm(0, 0, one(T), ShiftTerm(λ),
                                        Vector{OrbitalIntegral{<:Any,O,T,B,OV}}())
    OrbitalHamiltonian{O,T,B,OV,Proj}(h.R, vcat(h.terms, shift_term),
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
SCF.KrylovWrapper(hamiltonian::OrbitalHamiltonian{O,T,B,OV,Proj}) where {O,T,B,OV,Proj} =
    KrylovWrapper{T,OrbitalHamiltonian{O,T,B,OV,Proj}}(hamiltonian)

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
    copyto!(A.hamiltonian.R*y, A.hamiltonian⋆(A.hamiltonian.R*x))
