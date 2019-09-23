# * Orbital integrals

"""
    OrbitalIntegral{N}

Abstract type for integrals of rank `N` of orbitals, whose values need
to be recomputed every time the orbitals are updated. Rank 0
corresponds to a scalar value, rank 1 to a diagonal matrix, etc.
"""
abstract type OrbitalIntegral{N,aO<:AbstractOrbital,bO<:AbstractOrbital,T} end

Base.iszero(::OrbitalIntegral) = false

# ** Zero integral
struct ZeroIntegral{aO,bO,T} <: OrbitalIntegral{0,aO,bO,T} end

integral_value(::ZeroIntegral{aO,bO,T}) where {aO,bO,T} = zero(T)

SCF.update!(::ZeroIntegral, ::Atom; kwargs...) = nothing

# ** Orbital overlap integral

"""
    OrbitalOverlapIntegral(a, b, av, bv, value)

Represents the orbital overlap integral `⟨a|b⟩`, for orbitals `a` and
`b`, along with `view`s of their radial orbitals `av` and `bv` and the
current `value` of the integral.
"""
mutable struct OrbitalOverlapIntegral{aO,bO,T,OV} <: OrbitalIntegral{0,aO,bO,T}
    a::aO
    b::bO
    av::OV
    bv::OV
    value::T
end

function OrbitalOverlapIntegral(a::aO, b::bO, atom::Atom{T}) where {aO,bO,T}
    oo = OrbitalOverlapIntegral(a, b, view(atom, a), view(atom, b), zero(T))
    SCF.update!(oo)
    oo
end

Base.show(io::IO, oo::OrbitalOverlapIntegral) =
    write(io, "⟨$(oo.a)|$(oo.b)⟩")

"""
    SCF.update!(oo::OrbitalOverlapIntegral)

Update the value of the integral `oo`.
"""
function SCF.update!(oo::OrbitalOverlapIntegral; kwargs...)
    oo.value = materialize(applied(*, oo.av', oo.bv))[1]
end

"""
    SCF.update!(oo::OrbitalOverlapIntegral)

Update the value of the integral `oo` with respect to `atom`.
"""
function SCF.update!(oo::OrbitalOverlapIntegral, atom::Atom; kwargs...)
    oo.av = view(atom, oo.a)
    oo.bv = view(atom, oo.b)
    SCF.update!(oo; kwargs...)
end

integral_value(oo::OrbitalOverlapIntegral) = oo.value

# ** Operator matrix element

"""
    OperatorMatrixElement(a, b, Â, coeff, value)

Represents the matrix element `coeff*⟨a|Â|b⟩`, for the operator `Â`
and the orbitals `a` and `b`, along with the current `value` of the
integral. Typically, `Â` is a radial part of an operator and `coeff`
is the associated angular coefficient; `coeff` can be of any type
convertible to a scalar.
"""
mutable struct OperatorMatrixElement{aO,bO,OV,RO,T,Coeff} <: OrbitalIntegral{0,aO,bO,T}
    a::aO
    b::bO
    av::OV
    bv::OV
    Â::RO
    coeff::Coeff
    value::T
end

function OperatorMatrixElement(a::aO, b::bO, Â::RO, atom::Atom{T}, coeff) where {aO,bO,T,RO}
    ome = OperatorMatrixElement(a, b, view(atom, a), view(atom, b),
                                Â, coeff, zero(T))
    SCF.update!(ome)
    ome
end

function OperatorMatrixElement(a::aO, b::bO, A::M, atom::Atom{T}, coeff) where {aO,bO,T,
                                                                                M<:AbstractMatrix}
    R = radial_basis(atom)
    OperatorMatrixElement(a, b, applied(*,R,A,R'), atom, coeff)
end

function Base.show(io::IO, ome::OperatorMatrixElement)
    write(io, "⟨$(ome.a)|")
    show(io, ome.Â)
    write(io, "|$(ome.b)⟩")
end

function SCF.update!(ome::OperatorMatrixElement{aO,bO,OV,RO,T}; kwargs...) where {aO,bO,OV,RO,T}
    # NB: It is assumed that ome.Â, if necessary, is updated /before/
    # SCF.update!(ome) is called.
    #
    # TODO: This is not particularly efficient, since it allocates a
    # temporary vector; rewrite using applied(*, ...).
    b = ome.bv
    tmp = similar(b)
    tmp.args[2] .= zero(T)
    materialize!(MulAdd(1, ome.Â, b, 0, tmp))
    ome.value = convert(T,ome.coeff)*materialize(applied(*, ome.av', tmp))
end

function SCF.update!(ome::OperatorMatrixElement, atom::Atom; kwargs...)
    ome.av = view(atom, ome.a)
    ome.bv = view(atom, ome.b)
    SCF.update!(ome; kwargs...)
end

integral_value(ome::OperatorMatrixElement) = ome.value

# ** Hartree–Fock potentials

const HFPotentialOperator{T,B<:Basis} = RadialOperator{T,B,Diagonal{T,Vector{T}}}

"""
    HFPotential(k, a, b, av, bv, V̂, poisson)

Represents the `k`:th multipole exansion of the Hartree–Fock potential
formed by orbitals `a` and `b` (`av` and `bv` being `view`s of their
corresponding radial orbitals). `V̂` is the resultant one-body
potential formed, which can act on a third orbital and `poisson`
computes the potential by solving Poisson's problem.
"""
mutable struct HFPotential{kind,aO,bO,T,
                           OV<:RadialOrbital,
                           RO<:HFPotentialOperator{T},P<:AbstractPoissonProblem} <: OrbitalIntegral{1,aO,bO,T}
    k::Int
    a::aO
    b::bO
    av::OV
    bv::OV
    V̂::RO
    poisson::P
end
HFPotential(kind::Symbol, k::Int, a::aO, b::bO, av::OV, bv::OV, V̂::RO, poisson::P) where {aO,bO,T,OV,RO<:RadialOperator{T},P} =
    HFPotential{kind,aO,bO,T,OV,RO,P}(k, a, b, av, bv, V̂, poisson)

get_poisson(::CoulombInteraction{Nothing}, args...; kwargs...) =
    PoissonProblem(args...; kwargs...)

get_poisson(g::CoulombInteraction{<:AbstractQuasiMatrix}, args...; kwargs...) =
    AsymptoticPoissonProblem(args..., g.o; kwargs...)

function HFPotential(kind::Symbol, k::Int, a::aO, b::bO, atom::Atom{T}, g::CoulombInteraction; kwargs...) where {aO,bO,T}
    av, bv = view(atom, a), view(atom, b)
    R = av.args[1]
    D = Diagonal(Vector{T}(undef, size(R,2)))
    D.diag .= zero(T)
    V̂ = applied(*, R, D, R')
    poisson = get_poisson(g, k, av, bv; w′=applied(*, R, D.diag), kwargs...)
    update!(HFPotential(kind, k, a, b, av, bv, V̂, poisson), atom)
end

Base.convert(::Type{HFPotential{kind,aO₁,bO₁,T,OV,RO,P}},
             hfpotential::HFPotential{kind,aO₂,bO₂,T,OV,RO,P}) where {kind,aO₁,bO₁,aO₂,bO₂,T,OV,RO,P} =
                 HFPotential{kind,aO₁,bO₁,T,OV,RO,P}(hfpotential.k,
                                                     hfpotential.a, hfpotential.b,
                                                     hfpotential.av, hfpotential.bv,
                                                     hfpotential.V̂, hfpotential.poisson)

# *** Direct potential

"""
    DirectPotential

Special case of [`HFPotential`](@ref) for the direct interaction, in
which case the potential formed from two orbitals can be precomputed
before acting on a third orbital.
"""
const DirectPotential{aO,bO,T,OV,RO,P} = HFPotential{:direct,aO,bO,T,OV,RO,P}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.a), $(Y.b))")

"""
    SCF.update!(p::DirectPotential)

Update the direct potential `p` by solving the Poisson problem with
the current values of the orbitals forming the mutual density.
"""
function SCF.update!(p::DirectPotential{aO,bO,T,OV,RO,P}; kwargs...) where {aO,bO,T,OV,RO,P}
    p.poisson(p.av .⋆ p.bv; kwargs...)
    p
end

"""
    SCF.update!(p::DirectPotential, atom::Atom)

Update the direct potential `p` by solving the Poisson problem with
the current values of the orbitals of `atom` forming the mutual
density.
"""
function SCF.update!(p::DirectPotential{aO,bO,T,OV,RO,P}, atom::Atom; kwargs...) where {aO,bO,T,OV,RO,P}
    p.av = view(atom, p.a)
    p.bv = view(atom, p.b)
    SCF.update!(p; kwargs...)
end

"""
    materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest})

Materialize the lazy multiplication–addition of the type `y ←
α*V̂*x + β*y` where `V̂` is a [`DirectPotential`](@ref) (with a
precomputed direct potential computed via `SCF.update!`) and `x` and
`y` are [`RadialOrbital`](@ref)s.
"""
LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest}) where {T,Source,Dest} =
    materialize!(MulAdd(ma.α, ma.A.V̂, ma.B, ma.β, ma.C))

# *** Exchange potential

"""
    ExchangePotential

Special case of [`HFPotential`](@ref) for the exchange interaction, in
which case the potential is formed from the orbital acted upon, along
with another orbital, and then applied to a third orbital. Thus this
potential *cannot* be precomputed, but must be recomputed every time
the operator is applied. This makes this potential expensive to handle
and the number of times it is applied should be minimized, if possible.
"""
const ExchangePotential{aO,bO,T,OV,RO,P} = HFPotential{:exchange,aO,bO,T,OV,RO,P}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.b)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.a), ●)")

SCF.update!(p::ExchangePotential; kwargs...) = p

function SCF.update!(p::ExchangePotential, atom::Atom; kwargs...)
    p.av = view(atom, p.a)
    p.bv = view(atom, p.b)
    p
end

"""
    materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest})

Materialize the lazy multiplication–addition of the type `y ← α*V̂*x +
β*y` where `V̂` is a [`ExchangePotential`](@ref) (by solving the
Poisson problem with `x` as one of the constituent source orbitals in
the mutual density) and `x` and `y` are [`RadialOrbital`](@ref)s.
"""
function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest}) where {T,Source,Dest}
    p = ma.A
    p.poisson(p.av .⋆ ma.B) # Form exchange potential from conj(p.a)*b
    # Act with the exchange potential on p.bv
    materialize!(MulAdd(ma.α, p.V̂, p.bv, ma.β, ma.C))
end

# * Source terms

"""
    SourceTerm(operator, source_orbital, ov)

The point of `SourceTerm` is to implement inhomogeneous terms that
contribute to the equation for an orbital, and whose input is some
other `source_orbital`. This kind of term appears in
multi-configurational problems.
"""
mutable struct SourceTerm{QO,O,OV}
    operator::QO
    source_orbital::O
    ov::OV
end

Base.show(io::IO, st::SourceTerm) = write(io, "SourceTerm($(st.operator)|$(st.source_orbital)⟩)")

Base.iszero(::SourceTerm) = false
Base.similar(st::SourceTerm) = similar(st.ov)

function Base.copyto!(dest::Mul{<:Any,<:Tuple{<:AbstractQuasiMatrix,<:AbstractArray{<:Any,N}}},
                      src::Mul{<:Any,<:Tuple{<:AbstractQuasiMatrix,<:AbstractArray{<:Any,N}}}) where N
    d = last(dest.args)
    s = last(src.args)
    copyto!(IndexStyle(d), d, IndexStyle(s), s)
    dest
end

update!(st::SourceTerm; kwargs...) = nothing

function update!(st::SourceTerm, atom::Atom; kwargs...)
    st.ov = view(atom, st.source_orbital)
end

# These are source terms that do not depend on the atom (e.g. external
# source term, or a constant orbital).
update!(::SourceTerm{<:IdentityOperator,<:AbstractString}, ::Atom) = nothing

LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:SourceTerm, Source, Dest}) where {T,Source,Dest} =
    materialize!(MulAdd(ma.α, ma.A.operator, ma.A.ov, ma.β, ma.C))

function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:IdentityOperator{1}, Source, Dest}) where {T,Source,Dest}
    isone(ma.β) || lmul!(ma.β, ma.C.args[2])
    BLAS.axpy!(ma.α, ma.B.args[2], ma.C.args[2])
end

# * Shift terms

"""
    ShiftTerm(λ)

The point of `ShiftTerm` is to implement an overall energy shift of
the Hamiltonian.
"""
struct ShiftTerm{T}
    shift::UniformScaling{T}
end

Base.iszero(::ShiftTerm) = false

LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ShiftTerm, Source, Dest}) where {T,Source,Dest} =
    BLAS.axpy!(ma.α*ma.A.shift.λ, ma.B.args[2], ma.C.args[2])
