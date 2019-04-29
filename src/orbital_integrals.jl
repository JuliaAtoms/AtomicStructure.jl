# * Orbital integrals

"""
    OrbitalIntegral{N}

Abstract type for integrals of rank `N` of orbitals, whose values need
to be recomputed every time the orbitals are updated. Rank 0
corresponds to a scalar value, rank 1 to a diagonal matrix, etc.
"""
abstract type OrbitalIntegral{N,O<:AbstractOrbital,T,B<:Basis,OV<:RadialOrbital{T,B}} end

# ** Orbital overlap integral

"""
    OrbitalOverlapIntegral(a, b, av, bv, value)

Represents the orbital overlap integral `⟨a|b⟩`, for orbitals `a` and
`b`, along with `view`s of their radial orbitals `av` and `bv` and the
current `value` of the integral.
"""
mutable struct OrbitalOverlapIntegral{O,T,B,OV} <: OrbitalIntegral{0,O,T,B,OV}
    a::O
    b::O
    av::OV
    bv::OV
    value::T
end

function OrbitalOverlapIntegral(a::O, b::O, av::OV, bv::OV) where {O,T,B<:Basis,OV<:RadialOrbital{T,B}}
    oo = OrbitalOverlapIntegral(a, b, av, bv, zero(T))
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
    oo.value = oo.av'oo.bv
end

integral_value(oo::OrbitalOverlapIntegral) = oo.value

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
mutable struct HFPotential{kind,O,T,B<:Basis,OV,
                           RO<:HFPotentialOperator{T,B},P<:PoissonProblem} <: OrbitalIntegral{1,O,T,B,OV}
    k::Int
    a::O
    b::O
    av::OV
    bv::OV
    V̂::RO
    poisson::P
end
HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV, V̂::RO, poisson::P) where {O,T,B<:Basis,OV<:RadialOrbital{T,B},RO<:RadialOperator{T,B},P} =
    HFPotential{kind,O,T,B,OV,RO,P}(k, a, b, av, bv, V̂, poisson)

function HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV) where {O,T,B<:Basis,OV<:RadialOrbital{T,B}}
    R = av.args[1]
    D = Diagonal(Vector{T}(undef, size(R,2)))
    D.diag .= zero(T)
    V̂ = applied(*, R, D, R')
    poisson = PoissonProblem(k, av, bv, w′=applied(*, R, D.diag))
    update!(HFPotential(kind, k, a, b, av, bv, V̂, poisson))
end

Base.convert(::Type{HFPotential{kind,O₁,T,B,OV,RO}},
             hfpotential::HFPotential{kind,O₂,T,B,OV,RO,P}) where {kind,O₁,O₂,T,B,OV,RO,P} =
                 HFPotential{kind,O₁,T,B,OV,RO,P}(hfpotential.k,
                                                  hfpotential.a, hfpotential.b, hfpotential.av, hfpotential.bv,
                                                  hfpotential.V̂, hfpotential.poisson)

# *** Direct potential

"""
    DirectPotential

Special case of [`HFPotential`](@ref) for the direct interaction, in
which case the potential formed from two orbitals can be precomputed
before acting on a third orbital.
"""
const DirectPotential{O,T,B,OV,RO,P} = HFPotential{:direct,O,T,B,OV,RO,P}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.a), $(Y.b))")

"""
    SCF.update!(p::DirectPotential)

Update the direct potential `p` by solving the Poisson problem with
the current values of the orbitals forming the mutual density.
"""
function SCF.update!(p::DirectPotential{O,T,B,OV,RO,P}; kwargs...) where {O,T,B,OV,RO,P}
    p.poisson(;kwargs...)
    p
end

"""
    materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest})

Materialize the lazy multiplication–addition of the type `y ←
α*V̂*x + β*y` where `V̂` is a [`DirectPotential`](@ref) (with a
precomputed direct potential computed via `SCF.update!`) and `x` and
`y` are [`RadialOrbital`](@ref)s.
"""
LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest}) where {T,Source,Dest} =
    materialize!(MulAdd(ma.α, ma.A.V̂.args[2], ma.B.args[2],
                        ma.β, ma.C.args[2]))

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
const ExchangePotential{O,T,B,OV,RO,P} = HFPotential{:exchange,O,T,B,OV,RO,P}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.b)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.a), ●)")

# We can't update the exchange potentials, since they depend on the
# orbital they act on.
SCF.update!(p::ExchangePotential; kwargs...) = p

"""
    materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest})

Materialize the lazy multiplication–addition of the type `y ← α*V̂*x +
β*y` where `V̂` is a [`ExchangePotential`](@ref) (by solving the
Poisson problem with `x` as one of the constituent source orbitals in
the mutual density) and `x` and `y` are [`RadialOrbital`](@ref)s.
"""
function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest}) where {T,Source,Dest}
    p = ma.A
    p.poisson(ma.B) # Form exchange potential from conj(p.a)*b
    # Act with the exchange potential on p.bv
    materialize!(MulAdd(ma.α, p.V̂.args[2], p.bv.args[2],
                        ma.β, ma.C.args[2]))
end

# * Source terms

"""
    SourceTerm(operator, source_orbital, ov)

The point of `SourceTerm` is to implement inhomogeneous terms that
contribute to the equation for an orbital, and whose input is some
other `source_orbital`. This kind of term appears in
multi-configurational problems.
"""
struct SourceTerm{QO,O,OV}
    operator::QO
    source_orbital::O
    ov::OV
end

LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:SourceTerm, Source, Dest}) where {T,Source,Dest} =
    materialize!(MulAdd(ma.α, ma.A.operator, ma.A.ov, ma.β, ma.C))

function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:IdentityOperator{1}, Source, Dest}) where {T,Source,Dest}
    lmul!(ma.β, ma.C.args[2])
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

LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ShiftTerm, Source, Dest}) where {T,Source,Dest} =
    BLAS.axpy!(ma.α*ma.A.shift.λ, ma.B.args[2], ma.C.args[2])
