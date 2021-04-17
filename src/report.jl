function shift_unit(u::U, d) where {U<:Unitful.Unit}
    d = Int(3floor(Int, d/3))
    iszero(d) && return u
    for tens = Unitful.tens(u) .+ (d:(-3*sign(d)):0)
        haskey(Unitful.prefixdict, tens) && return U(tens, u.power)
    end
    u
end

function shift_unit(u::Unitful.FreeUnits, d)
    tu = typeof(u)
    us,ds,a = tu.parameters

    uu = shift_unit(us[1], d)

    Unitful.FreeUnits{(uu,us[2:end]...), ds, a}()
end

function si_round(q::Quantity; fspec="{1:+9.4f} {2:s}")
    v,u = ustrip(q), unit(q)
    if !iszero(v)
        u = shift_unit(u, log10(abs(v)))
        q = u(q)
    end
    format(fspec, ustrip(q), unit(q))
end

function report(io::IO, fock::Fock{<:Atom{T}}; powers = (2,1,-1,-2,-3)) where T
    atom = fock.quantum_system
    Z = charge(atom.potential)
    N = num_electrons(atom)
    println(io, "Atom{$(T)}")
    println(io, "R = ", radial_basis(atom))
    show(io, atom.potential)
    println(io, ", $(N) e⁻ ⇒ Q = $(Z-N)")

    Eₜₒₜ = SCF.energy(fock, :total_energy)*u"hartree"
    Eₖᵢₙ = SCF.energy(fock, :kinetic_energy)*u"hartree"
    Eₚₒₜ = Eₜₒₜ - Eₖᵢₙ
    VT = Eₚₒₜ/Eₖᵢₙ

    # At some point we will support comparing to reference data, taken from e.g.
    # - Saito, S. L. (2009). Hartree-Fock-Roothaan energies and expectation
    #   values for the neutral atoms he to uuo: the B-spline expansion
    #   method. Atomic Data and Nuclear Data Tables, 95(6),
    #   836–870. http://dx.doi.org/10.1016/j.adt.2009.06.001
    total_energy_errors = Any["", "", "", VT+2]

    pretty_table(io, hcat(["Total energy", "Kinetic energy", "Potential energy", "Virial ratio"],
                          vcat(si_round.((Eₜₒₜ, Eₖᵢₙ, Eₚₒₜ), fspec="{1:+15.10f} {2:s}")..., format("{1:+15.10f}", VT)),
                          total_energy_errors),
                 noheader=true,
                 alignment=:l,tf=tf_borderless)

    ϵ = si_round.((SCF.orbital_energies(fock) ./ degeneracy.(atom.orbitals)) * u"hartree")

    R = radial_basis(atom)
    r = axes(R,1)

    rᵏ = collect(R'*QuasiDiagonal(r.^k)*R for k in powers)
    radial_moments = zeros(T, length(atom.orbitals), length(powers))
    for i = eachindex(atom.orbitals)
        ϕ = view(atom.radial_orbitals.args[2], :, i)
        for (j,rᵏ) in enumerate(rᵏ)
            radial_moments[i,j] = dot(ϕ, rᵏ, ϕ)
        end
    end

    # Should also compute electron density ρ(0) and Kato cusp

    println(io)

    pretty_table(io, [atom.orbitals ϵ radial_moments],
                 header=vcat(["Orbital", "ϵ"], ["⟨r"*(k≠1 ? to_superscript(k) : "")*"⟩" for k in powers]),
                 tf=tf_borderless)
end

report(fock::Fock; kwargs...) = report(stdout, fock; kwargs...)
