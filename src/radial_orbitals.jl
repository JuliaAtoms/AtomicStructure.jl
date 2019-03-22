"""
    RadialOrbital

A radial orbital is represented using basis coupled with a vector of
expansion coefficients with respect to that basis; the basis
implemented as an `AbstractQuasimatrix`.
"""
const RadialOrbital{T,B<:AbstractQuasiMatrix} = Mul{<:Any,<:Tuple{<:B,<:AbstractVector{T}}}
const AdjointRadialOrbital{T,B<:AbstractQuasiMatrix} = Mul{<:Any,<:Tuple{<:Adjoint{<:Any,<:AbstractVector{T}},<:QuasiAdjoint{<:Any,<:B}}}

"""
    RadialOrbitals

A collection of radial orbitals is instead represented using a
matrix of such expansion coefficients, where each matrix column
corresponds to a single radial orbital.
"""
const RadialOrbitals{T,B<:AbstractQuasiMatrix} = Mul{<:Any,<:Tuple{<:B,<:AbstractMatrix{T}}}

const RadialOperator{T,B<:AbstractQuasiMatrix,M<:AbstractMatrix{T}} =
    Mul{<:Any,<:Tuple{<:B,M,<:QuasiAdjoint{<:Any,<:B}}}

matrix(o::RadialOperator) = o.args[2]

function Base.zero(o::RadialOperator)
    A,B,C = o.args
    applied(*, A, Diagonal(zeros(eltype(B),size(B,1))), C)
end

function Base.:(+)(a::RadialOperator{T,B}, b::RadialOperator{T,B}) where {T,B}
    R = first(a.args)
    applied(*, R, matrix(a) + matrix(b), R')
end
