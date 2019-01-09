"""
    hf_equations(csf[; verbosity=0])

Derive the Hartree–Fock equations for the non-relativistic
configuration state function `csf`. All constituent orbitals are
assumed spectroscopic, i.e. they are all assigned Lagrange multipliers
to ensure their orthonormality. Additionally, orbitals of the same `ℓ`
but different `n` are assigned off-diagonal Lagrange multipliers.
"""
function hf_equations(csf::NonRelativisticCSF; verbosity=0)
    verbosity > 0 &&
        println("Deriving HF equations for $(csf)")
    
    orbitals = csf.config.orbitals
    eng = energy_expression(csf; verbosity=verbosity)[2][1]
    λs = lagrange_multipliers(orbitals)
    eng += λs
    
    if verbosity > 2
        println("Orbitals: $(orbitals)")
        println("Lagrange multipliers: $(λs)")
    end
    verbosity > 1 && println("E = $eng\n")

    map(peel(csf.config)) do(orb,occ,state)
        hf_eqs = diff(eng, orb, occ)
        verbosity > 2 && println("∂[$(orb)] E = $(hf_eqs)")
        hf_eqs
    end
end

function Base.diff(atom::Atom; kwargs...)
    length(atom.csfs) > 1 &&
        throw(ArgumentError("Cannot derive Hartree–Fock equations for a multi-configurational atom"))
    
    hf_equations(atom.csfs[1]; kwargs...)
end

Base.view(fock::Fock{A,E}, args...) where {A<:Atom,E} =
    view(fock.quantum_system.radial_orbitals.mul.factors[2], args...)

function LinearAlgebra.ldiv!(fock::Fock{A,E}, c::M) where {A<:Atom,E,M<:AbstractMatrix}
    c .= rand() # Bogus solution of secular equations
end
