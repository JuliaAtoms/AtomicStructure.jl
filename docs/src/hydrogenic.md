# Hydrogenic initialization

There are two modes of hydrogenic initialization implemented, plain
hydrogenics, where all orbitals are initialized to as eigenorbitals of
the bare charge ``Z``, and screened hydrogenics, where each orbital
``i`` is initialized as the eigenorbital of an efficient charge
``Z-\sigma_i``, where ``\sigma_i`` to some extent accounts for the
screening of orbital ``i`` due to all the other electrons in the
configuration.

```@meta
CurrentModule = AtomicStructure
DocTestSetup = quote
    using AtomicStructure
    using AtomicLevels
end
```

## Plain hydrogenics

```@docs
hydrogenic!
```

## Screened hydrogenics

```@docs
screened_hydrogenic!
AtomicStructure.screening
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```
