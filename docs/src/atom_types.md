# Atom types

```@meta
CurrentModule = AtomicStructure
DocTestSetup = quote
    using AtomicStructure
end
```

```@docs
Atom
Atom(::UndefInitializer, ::Type{T}, R::B, configurations::Vector{TC}, potential::P, ::Type{C}, mix_coeffs::CV=vcat(one(C), zeros(C, length(configurations)-1))) where {T<:Number,B<:BasisOrRestricted,TC,C,CV<:AbstractVector{<:C},P}
Atom(init::Symbol, ::Type{T}, R::B, args...; kwargs...) where {T<:Number,B<:BasisOrRestricted,TC,C,P}
Atom(init::Init, R::B, args...; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T}}
Atom(R::B, configurations::Vector{<:TC}, potential::P, ::Type{C}=eltype(R), args...; kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,C,P<:AbstractPotential}
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
CurrentModule = AtomicStructure
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
