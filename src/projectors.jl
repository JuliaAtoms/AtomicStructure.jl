"""
    Projector(ϕs, orbitals)

Represents the projector on the subspace spanned by the radial
orbitals `ϕs` (corresponding to `orbitals`).
"""
struct Projector{T,B<:AbstractQuasiMatrix,RO<:RadialOrbital{T,B},O} <: NBodyOperator{1}
    ϕs::Vector{RO}
    orbitals::Vector{O}
end

Base.iszero(me::OrbitalMatrixElement{1,A,<:Projector,B}) where {A<:SpinOrbital,B<:SpinOrbital} =
    me.a[1] ∉ me.o.orbitals || me.b[1] ∉ me.o.orbitals

function Base.show(io::IO, projector::Projector)
    write(io, "P(")
    write(io, join(string.(projector.orbitals), " "))
    write(io, ")")
end

projectout!(y::RO, ::Nothing) where RO = y

"""
    projectout!(y, projector)

Project out all components of `y` parallel to the radial orbitals
`projector.ϕs`.
"""
function projectout!(y::RO, projector::Proj) where {RO,Proj<:Projector}
    yc = y.args[2]
    
    for ϕ in projector.ϕs
        c = materialize(applied(*, ϕ', y))
        yc .-= c*ϕ.args[2]
        # y -= c*ϕ
    end
    y
end
