# Equation systems

```@meta
CurrentModule = Atoms
DocTestSetup = quote
    using Atoms
end
```

```@docs
AtomicEquations
SCF.update!(equations::AtomicEquations; kwargs...)
SCF.energy_matrix!(H::HM, hfeqs::AtomicEquations, which::Symbol=:total) where {HM<:AbstractMatrix}
find_symmetries
generate_atomic_orbital_equations
Base.diff
```

## Common integrals

When deriving the equations of motion from an energy expression, the
same integral may appear many times in the equations for different
orbitals, multiplied by different factors, etc. To minimize the
reevaluation of integrals, [`AtomicEquations`](@ref) keeps track of
all the common integrals, and they are recomputed exactly once, when
[`SCF.update!`](@ref) is called. The routines below are used when
setting up the equation system.

```@docs
pushterms!
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
