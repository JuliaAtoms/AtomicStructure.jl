# * Orbital integrals

"""
    OrbitalIntegral{N}

Abstract type for integrals of rank `N` of orbitals, whose values need
to be recomputed every time the orbitals are updated. Rank 0
corresponds to a scalar value, rank 1 to a diagonal matrix, etc.
"""
abstract type OrbitalIntegral{N,O<:AbstractOrbital,T,B<:AbstractQuasiMatrix,OV<:RadialOrbital{T,B}} end

# ** Orbital overlap integral

mutable struct OrbitalOverlapIntegral{O,T,B,OV} <: OrbitalIntegral{0,O,T,B,OV}
    a::O
    b::O
    av::OV
    bv::OV
    value::T
end

function OrbitalOverlapIntegral(a::O, b::O, av::OV, bv::OV) where {O,T,B,OV}
    oo = OrbitalOverlapIntegral(a, b, av, bv, zero(T))
    SCF.update!(oo)
    oo
end

Base.show(io::IO, oo::OrbitalOverlapIntegral) =
    write(io, "⟨$(oo.a)|$(oo.b)⟩")

function SCF.update!(oo::OrbitalOverlapIntegral; kwargs...)
    oo.value = oo.av'oo.bv
end

integral_value(oo::OrbitalOverlapIntegral) = oo.value

# ** Hartree–Fock potentials

const HFPotentialOperator{T,B} = RadialOperator{T,B,Diagonal{T,Vector{T}}}

mutable struct HFPotential{kind,O,T,B,OV,
                           RO<:HFPotentialOperator{T,B},P<:PoissonProblem} <: OrbitalIntegral{1,O,T,B,OV}
    k::Int
    a::O
    b::O
    av::OV
    bv::OV
    V̂::RO
    poisson::P
end
HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV, V̂::RO, poisson::P) where {O,T,B,OV<:RadialOrbital{T,B},RO<:RadialOperator{T,B},P} =
    HFPotential{kind,O,T,B,OV,RO,P}(k, a, b, av, bv, V̂, poisson)

function HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV) where {O,T,B,OV<:RadialOrbital{T,B}}
    R = av.mul.factors[1]
    D = Diagonal(Vector{T}(undef, size(R,2)))
    D.diag .= zero(T)
    V̂ = R*D*R'
    poisson = PoissonProblem(k, av, bv, w′=R*D.diag)
    update!(HFPotential(kind, k, a, b, av, bv, V̂, poisson))
end

Base.convert(::Type{HFPotential{kind,O₁,T,B,OV,RO}},
             hfpotential::HFPotential{kind,O₂,T,B,OV,RO,P}) where {kind,O₁,O₂,T,B,OV,RO,P} =
                 HFPotential{kind,O₁,T,B,OV,RO,P}(hfpotential.k,
                                                  hfpotential.a, hfpotential.b, hfpotential.av, hfpotential.bv,
                                                  hfpotential.V̂, hfpotential.poisson)

# *** Direct potential

const DirectPotential{O,T,B,OV,RO,P} = HFPotential{:direct,O,T,B,OV,RO,P}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.a), $(Y.b))")

function SCF.update!(p::DirectPotential{O,T,B,OV,RO,P}; kwargs...) where {O,T,B,OV,RO,P}
    p.poisson(;kwargs...)
    p
end

LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest}) where {T,Source,Dest} =
    materialize!(MulAdd(ma.α, ma.A.V̂.mul.factors[2], ma.B.mul.factors[2],
                        ma.β, ma.C.mul.factors[2]))

# *** Exchange potential

const ExchangePotential{O,T,B,OV,RO,P} = HFPotential{:exchange,O,T,B,OV,RO,P}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.b)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.a), ●)")

# We can't update the exchange potentials, since they depend on the
# orbital they act on.
SCF.update!(p::ExchangePotential; kwargs...) = p

function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest}) where {T,Source,Dest}
    p = ma.A
    p.poisson(ma.B) # Form exchange potential from conj(p.a)*b
    # Act with the exchange potential on p.bv
    materialize!(MulAdd(ma.α, p.V̂.mul.factors[2], p.bv.mul.factors[2],
                        ma.β, ma.C.mul.factors[2]))
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
