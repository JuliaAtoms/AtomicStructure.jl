import AtomicLevels: spectroscopic_label

"""
    hydrogenic!(atom[; kwargs...])

Initialize the radial orbitals of `atom` to their hydrogenic
values. This is done via simple diagonalization of the one-body
Hamiltonian for each angular symmetry. The `kwargs` are passed on to
`diagonalize_one_body` and can be used to influence how the
diagonalization is performed.
"""
function hydrogenic!(atom::Atom{T,T,B,O,TC,C,CM,P}; verbosity=0, kwargs...) where {T,B,O,TC,C,CM,P}
    # This is only for non-relativistic atoms
    all_bound(atom) ||
        throw(ArgumentError("Cannot initialize $(atom) to hydrogenic state if not all orbitals are bound"))

    print_block() do io
        verbosity > 0 && println(io, "Hydrogenic initialization of the orbitals of $(atom)")

        R,Φ = atom.radial_orbitals.mul.factors
        Φ .= zero(T)

        # Orbitals are grouped by ℓ, since they have the same one-body
        # Hamiltonian. They should rather be grouped by
        # `symmetry(orb)`, where `symmetry` dispatches on the type of
        # the orbital. This would allow for pseudo-potentials that are
        # ℓ- and s-dependent.
        orbitals = Dict{Int,Vector{O}}()
        for orb in atom.orbitals
            orbitals[orb.ℓ] = vcat(get(orbitals, orb.ℓ, O[]), orb)
        end

        Z = charge(atom.potential)
        hydrogen_like_eng = n -> -Z^2/(2n^2)

        eng_fmt = FormatSpec("+10.7f")
        Δeng_fmt = FormatSpec("+10.3e")
        fmt_eng_vec(fspec, v) = "["*join(fmt.(Ref(fspec), v), ", ")*"]"

        for ℓ in keys(orbitals)
            print_block(io) do io
                # For shift-and-invert type of diagonalization, we target an
                # eigenvalue slightly below the lowest energy in partial wave
                # ℓ of a hydrogen-like system with point charge Z.
                min_n = ℓ + 1
                Iₚℓ = hydrogen_like_eng(min_n)
                σ = 2Iₚℓ

                max_n = maximum(o -> o.n, orbitals[ℓ])
                nev = max_n - ℓ

                verbosity > 1 &&
                    println(io, "Diagonalizing symmetry ℓ = $(spectroscopic_label(ℓ)), maximum n = $(max_n) => $(nev) eigenvalues required")
                verbosity > 2 &&
                    println(io, "Target eigenvalue: ≤ $(Iₚℓ) Ha")

                H = one_body_hamiltonian(atom, first(orbitals[ℓ]))
                λᴴ,Φᴴ = diagonalize_one_body(H, nev; σ=σ, io=io, verbosity=verbosity, kwargs...)
                for orb in orbitals[ℓ]
                    j = findfirst(isequal(orb), atom.orbitals)
                    copyto!(view(Φ, :, j), view(Φᴴ, :, orb.n-ℓ))
                end

                if verbosity > 2
                    λᴴₐ = hydrogen_like_eng.(min_n:max_n)
                    println(io, "Hydrogenic energies $(fmt_eng_vec(eng_fmt, λᴴ)) Ha")
                    println(io, "Analytic energies   $(fmt_eng_vec(eng_fmt, λᴴₐ)) Ha")
                    println(io, "Δ                   $(fmt_eng_vec(Δeng_fmt, λᴴ-λᴴₐ)) Ha")
                end
            end
        end

        print_block(io) do io
            for j in eachindex(atom.orbitals)
                # We index with the range j:j instead of the value j, since
                # the latter produces a ContinuumArrays.QuasiArrays.SubQuasiArray
                # which is tricky to dispatch on.
                N₀ = norm(atom.radial_orbitals[:,j:j])
                normalize!(atom.radial_orbitals[:,j:j])
                if verbosity > 2
                    N = norm(atom.radial_orbitals[:,j:j])
                    printfmtln(io, "Initial norm of {1:s}: {2:f}, 1-normalized: {3:e}",
                               atom.orbitals[j], N₀, 1-N)
                end
            end
        end
    end

    atom
end

export hydrogenic!