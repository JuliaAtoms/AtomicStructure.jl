const RestrictionMatrix = BandedMatrix{<:Int, <:FillArrays.Ones}

const RestrictedBasis{B<:Basis} = Mul{<:Any,<:Tuple{B, <:RestrictionMatrix}}
const AdjointRestrictedBasis{B<:Basis} = Mul{<:Any,<:Tuple{<:Adjoint{<:Any,<:RestrictionMatrix}, <:QuasiAdjoint{<:Any,B}}}

const RestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:RestrictedBasis{B}}
const AdjointRestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:AdjointRestrictedBasis{B}}

const BasisOrRestricted{B<:Basis} = Union{B,RestrictedBasis{<:B},<:RestrictedQuasiArray{<:Any,<:Any,<:B}}
const AdjointBasisOrRestricted{B<:Basis} = Union{<:QuasiAdjoint{<:Any,B},AdjointRestrictedBasis{<:B},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:B}}
