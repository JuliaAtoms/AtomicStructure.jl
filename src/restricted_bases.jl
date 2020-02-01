const RestrictionMatrix = BandedMatrix{<:Int, <:FillArrays.Ones}
const RestrictedQuasiArray{T,N,B<:Basis} = SubQuasiArray{T,N,B}
const AdjointRestrictedQuasiArray{T,N,B<:Basis} = QuasiAdjoint{T,<:RestrictedQuasiArray{T,N,B}}

const BasisOrRestricted{B<:Basis} = Union{B,<:RestrictedQuasiArray{<:Any,<:Any,<:B}}
const AdjointBasisOrRestricted{B<:Basis} = Union{<:QuasiAdjoint{<:Any,B},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:B}}
