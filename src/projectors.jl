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

function projectout!(y::AbstractVector{T}, projector::Projector) where T
    s = zero(T)
    for ϕ in projector.ϕs
        c = dot(ϕ.args[2], projector.S, y)
        y .-= c .* ϕ.args[2]
        s += c
    end
    s
end

function projectout!(Y::AbstractMatrix{T}, projector::Projector) where T
    s = zero(T)
    n = size(Y, 2)
    for j = 1:n
        s += projectout!(view(Y, :, j), projector)
    end
    s
end

"""
    projectout!(y, projector)

Project out all components of `y` parallel to the radial orbitals
`projector.ϕs`.
"""
projectout!(y::RadialOrbital, projector::Projector) =
    projectout!(y.args[2], projector)

function LinearAlgebra.Matrix(p::Projector{T}) where T
    n = size(p.S,1)
    U = zeros(T, n, length(p.ϕs))
    for (j,ϕ) in enumerate(p.ϕs)
        U[:,j] = ϕ.args[2]
    end
    adjU = U'p.S
    U*adjU
end

function LinearAlgebra.mul!(y::AbstractVector, p::Projector, x::AbstractVector,
                            α::Number=true, β::Number=false)
    if iszero(β)
        y .= false
    elseif !isone(β)
        lmul!(β, y)
    end
    iszero(α) && return y

    for ϕ in p.ϕs
        ϕv = ϕ.args[2]
        c = α*dot(ϕv, p.S, x)
        y .+= c .* ϕv
    end

    y
end

function LinearAlgebra.mul!(Y::AbstractMatrix, p::Projector, X::AbstractMatrix,
                            α::Number=true, β::Number=false)
    n = size(Y,2)
    for j = 1:n
        mul!(view(Y, :, j), p, view(X, :, j), α, β)
    end
    Y
end
