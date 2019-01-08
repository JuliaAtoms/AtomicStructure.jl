mutable struct TwoComponent{T}
    P::T
    Q::T
end

Base.promote_type(::Type{T}, ::Type{TwoComponent{T}}) where T = TwoComponent{T}

Base.one(::Type{TwoComponent{T}}) where T = TwoComponent(one(T), one(T))
Base.zero(::Type{TwoComponent{T}}) where T = TwoComponent(zero(T), zero(T))
