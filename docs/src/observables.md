# Observables

```@meta
CurrentModule = Atoms
DocTestSetup = quote
    using Atoms
end
```

```@docs
Observable
Observable(operator::QuantumOperator, atom::A, overlaps::Vector{<:OrbitalOverlap}, integrals::Vector{OrbitalIntegral}, integral_map::Dict{Any,Int}, symmetries::Dict) where {A<:Atom}
observe!
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
