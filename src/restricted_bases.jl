const RestrictionMatrix = BandedMatrix{<:Int, <:FillArrays.Ones}

const RestrictionTuple{B<:Basis} = Tuple{B, <:RestrictionMatrix}
const AdjointRestrictionTuple{B<:Basis} = Tuple{<:Adjoint{<:Any,<:RestrictionMatrix}, <:QuasiAdjoint{<:Any,B}}

const RestrictedBasis{B<:Basis} = Mul{<:Any,<:RestrictionTuple{B}}
const AdjointRestrictedBasis{B<:Basis} = Mul{<:Any,<:AdjointRestrictionTuple{B}}

const RestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:RestrictionTuple{B}}
const AdjointRestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:AdjointRestrictionTuple{B}}

const BasisOrRestricted{B<:Basis} = Union{B,RestrictedBasis{<:B},<:RestrictedQuasiArray{<:Any,<:Any,<:B}}
const AdjointBasisOrRestricted{B<:Basis} = Union{<:QuasiAdjoint{<:Any,B},AdjointRestrictedBasis{<:B},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:B}}
