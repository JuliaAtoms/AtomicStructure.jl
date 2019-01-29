# A radial orbital is represented using basis coupled with a vector of
# expansion coefficients with respect to that basis; the basis
# implemented as an `AbstractQuasimatrix`. A collection of radial
# orbitals is instead represented using a matrix of such expansion
# coefficients, where each matrix column corresponds to a single
# radial orbital.

const RadialCoeff{T<:Number} = Union{T,TwoComponent{T}}
const RadialOrbital{T<:RadialCoeff,B<:AbstractQuasiMatrix} = MulQuasiArray{T,1,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractVector}}}
const RadialOrbitals{T<:RadialCoeff,B<:AbstractQuasiMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,<:AbstractMatrix}}}

const RadialOperator{T<:RadialCoeff,B<:AbstractQuasiMatrix,M<:AbstractMatrix} = MulQuasiArray{T,2,<:Mul{<:Tuple,<:Tuple{<:B,M,<:QuasiAdjoint{<:Any,<:B}}}}
matrix(o::RadialOperator) = o.mul.factors[2]
function Base.zero(o::RadialOperator)
    A,B,C = o.mul.factors
    A*Diagonal(zeros(eltype(B),size(B,1)))*C
end

function Base.:(+)(a::RadialOperator{T,B}, b::RadialOperator{T,B}) where {T,B}
    R = first(a.mul.factors)
    R*(matrix(a) + matrix(b))*R'
end
