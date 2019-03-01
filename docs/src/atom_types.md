# Atom types

```@meta
CurrentModule = Atoms
DocTestSetup = quote
    using Atoms
end
```

```@docs
Atom
Atom(radial_orbitals::RadialOrbitals{T,B}, orbitals::Vector{O}, configurations::Vector{<:TC}, potential::P, ::Type{C}) where {T<:Number,B,O,TC<:ManyElectronWavefunction,C,P}
Atom(::UndefInitializer, ::Type{T}, R::B, configurations::Vector{TC}, potential::P, ::Type{C}) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
Atom(init::Symbol, ::Type{T}, R::B, configurations::Vector{TC}, potential::P, ::Type{C}; kwargs...) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}
Atom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C}; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,C,P<:AbstractPotential}
Atom(R::B, configurations::Vector{<:TC}, potential::P, ::Type{C}=eltype(R); kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,C,P<:AbstractPotential}
Atom(other_atom::Atom{T,B,O,TC,C,P}, configurations::Vector{<:TC}; kwargs...) where {T,B,O,TC,C,P}
DiracAtom
```

```@docs
getindex
view
num_electrons
```

## Internals

```@meta
CurrentModule = Atoms
```

```@docs
ManyElectronWavefunction
outsidecoremodel
all_bound
SCF.coefficients
SCF.orbitals
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
