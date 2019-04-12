"""
    Projector(ϕs)

Represents the projector *out of* the subspace spanned by the radial
orbitals `ϕs`
"""
struct Projector{T,B<:AbstractQuasiMatrix,RO<:RadialOrbital{T,B}}
    ϕs::Vector{RO}
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