# * Diagonal matrices

function diagmuladd!(α, A::Diagonal, B::AbstractVector, β, C::AbstractVector)
    mA, nA = size(A)
    mB = length(B)
    nA == mB || throw(DimensionMismatch("Dimensions must match"))
    length(C) == mA || throw(DimensionMismatch("Dimensions must match"))

    isone(β) || lmul!(β, C)
    (nA == 0 || mB == 0) && return C

    @inbounds for k = 1:mB
        C[k] += α * A.diag[k] * B[k]
    end

    C
end

function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T,
                                            <:RadialOperator{<:Any, <:Any, <:Diagonal},
                                            Source, Dest}) where {T,Source<:RadialOrbital,Dest<:RadialOrbital}
    diagmuladd!(ma.α, ma.A.args[2], ma.B.args[2], ma.β, ma.C.args[2])
    ma.C
end

ArrayLayouts.default_blasmul!(α, A::Diagonal, B::AbstractVector, β, C::AbstractVector) =
    diagmuladd!(α, A, B, β, C)

# * SymTridiagonal matrices

function symtridiagmuladd!(α, A::SymTridiagonal, B::AbstractVector, β, C::AbstractVector)
    mA, nA = size(A)
    m = length(B)
    mA == nA == m || throw(DimensionMismatch("Dimensions must match"))
    length(C) == m || throw(DimensionMismatch("Dimensions must match"))

    isone(β) || lmul!(β, C)
    (nA == 0 || m == 0) && return C

    d = A.dv
    e = A.ev

    @inbounds begin
        x₊ = B[1]
        x₀ = zero(x₊)
        # If m == 1 then e[1] is out of bounds
        e₀ = m > 1 ? zero(e[1]) : zero(eltype(e))
        for i = 1:m - 1
            x₋, x₀, x₊ = x₀, x₊, B[i + 1]
            e₋, e₀ = e₀, e[i]
            C[i] += α * (e₋*x₋ + d[i]*x₀ + e₀*x₊)
        end
        C[m] += α * (e₀*x₀ + d[m]*x₊)
    end

    C
end

function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T,
                                            <:RadialOperator{<:Any, <:Any, <:SymTridiagonal},
                                            Source, Dest}) where {T,Source<:RadialOrbital,Dest<:RadialOrbital}
    symtridiagmuladd!(ma.α, ma.A.args[2], ma.B.args[2], ma.β, ma.C.args[2])
    ma.C
end

ArrayLayouts.default_blasmul!(α, A::SymTridiagonal, B::AbstractVector, β, C::AbstractVector) =
    symtridiagmuladd!(α, A, B, β, C)
