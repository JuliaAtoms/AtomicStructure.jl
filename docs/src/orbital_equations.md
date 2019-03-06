# Orbital equations

```@meta
CurrentModule = Atoms
DocTestSetup = quote
    using Atoms
end
```

## Hamiltonian

```@docs
OrbitalHamiltonianTerm
coefficient
OrbitalHamiltonian
Projector
projectout!
SCF.energy_matrix!(H::HM, hamiltonian::OrbitalHamiltonian{O,T,B}, ϕ::RadialOrbital{T,B}) where {HM<:AbstractMatrix,O,T,B}
Base.filter(fun::Function, H::OrbitalHamiltonian)
Base.:(+)(h::OrbitalHamiltonian{O,T,B,OV,Proj}, λ::UniformScaling) where {O,T,B,OV,Proj}
Base.:(-)(h::OrbitalHamiltonian, λ::UniformScaling)
SCF.KrylovWrapper
```

## Orbital integrals and terms

```@docs
OrbitalIntegral
OrbitalOverlapIntegral
SCF.update!(oo::OrbitalOverlapIntegral; kwargs...)
HFPotential
DirectPotential
SCF.update!(p::DirectPotential{O,T,B,OV,RO,P}; kwargs...) where {O,T,B,OV,RO,P}
LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest}) where {T,Source,Dest}
ExchangePotential
LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest}) where {T,Source,Dest}
SourceTerm
ShiftTerm
```

## Orbital equations

```@docs
AtomicOrbitalEquation
energy(hfeq::AtomicOrbitalEquation, term=:all)
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
