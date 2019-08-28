# Observables

```@meta
CurrentModule = Atoms
DocTestSetup = quote
    using Atoms
end
```

```@docs
Observable
Observable(operator::QuantumOperator, atom::A, overlaps::Vector{<:OrbitalOverlap}, integrals::Vector{OrbitalIntegral}, integral_map::Dict{Any,Int}, symmetries::Dict, selector::Function; double_counted::Bool) where {T,A<:Atom{T}}
observe!
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
