"""
    Projector(ϕs, orbitals, S)

Represents the projector on the subspace spanned by the radial
orbitals `ϕs` (corresponding to `orbitals`).
"""
struct Projector{T,B<:AbstractQuasiMatrix,RO<:RadialOrbital{T,B},O,Metric} <: NBodyOperator{1}
    ϕs::Vector{RO}
    orbitals::Vector{O}
    S::Metric
end

Base.eltype(projector::Projector{T}) where T = T

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
function projectout!(y::RO, projector::Proj) where {RO,T,Proj<:Projector{T}}
    yc = y.args[2]

    s = zero(T)
    for ϕ in projector.ϕs
        c = dot(ϕ.args[2], projector.S, y.args[2])
        yc .-= c*ϕ.args[2]
        # y -= c*ϕ
        s += c
    end
    s
end

function LinearAlgebra.Matrix(p::Projector{T}) where T
    n = size(p.S,1)
    U = zeros(T, n, length(p.ϕs))
    for (j,ϕ) in enumerate(p.ϕs)
        U[:,j] = ϕ.args[2]
    end
    adjU = U'p.S
    U*adjU
end
