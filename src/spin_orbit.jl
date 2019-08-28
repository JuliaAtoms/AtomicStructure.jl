struct SpinOrbitInteraction <: OneBodyOperator end

atomic_hamiltonian(::Atom{T,B,O,TC,CV,P}) where {T,B,O,TC,CV,P<:Union{RelativisticPseudoPotential,PseudoPotentials.ScalarSORelativisticPseudoPotential}} =
    FieldFreeOneBodyHamiltonian() + CoulombInteraction() + SpinOrbitInteraction()

Base.show(io::IO, ::SpinOrbitInteraction) = write(io, "V̂ₛₒ")

"""
    mⱼ(o::SpinOrbital)

Return the total mⱼ that `o` couples to in LS coupling, i.e. `mℓ+mₛ`.
"""
mⱼ(o::SpinOrbital) =
    o.mℓ + (o.spin ? half(1) : half(-1))

function Base.iszero(me::OrbitalMatrixElement{1,A,SpinOrbitInteraction,B}) where {A<:SpinOrbital,B<:SpinOrbital}
    a = me.a[1]
    b = me.b[1]

    ℓa,ma = jmⱼ(a)
    ℓb,mb = jmⱼ(b)
    # a.orb.n ≠ b.orb.n ||
    ℓa ≠ ℓb || iszero(ℓa) || mⱼ(a) != mⱼ(b)
end


function get_operator(::SpinOrbitInteraction, atom::Atom, # {T,B,O,TC,CV,P},
                      a::aO, b::bO) where {aO,bO} # ,T,B,O,TC,CV,P<:RelativisticPseudoPotential}
    if a == b && atom.potential isa RelativisticPseudoPotential
        # In case of two-component relativistic pseudopotentials, the
        # orbital-diagonal spin–orbit contribution is already
        # accounted for.
        return 0
    end
    R = radial_basis(atom)
    r = locs(R)
    Vₛₒ = spin_orbit_potential(atom.potential, r, a, b)
    applied(*,R,Diagonal(Vₛₒ),R')
end

# get_operator(::SpinOrbitInteraction, ::Atom{T,B,O,TC,CV,P},
#              ::aO, ::bO) where {aO,bO,T,B,O,TC,CV,P<:AbstractPotential} =
#                  throw(ArgumentError("Spin–orbit interaction not yet implemented without RelativisticPseudoPotential"))
