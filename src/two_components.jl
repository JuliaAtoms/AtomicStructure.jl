mutable struct TwoComponent{T}
    P::T
    Q::T
end

Base.promote_type(::Type{T}, ::Type{TwoComponent{T}}) where T = TwoComponent{T}
