function LazyArrays.default_blasmul!(α, A::Diagonal, B::AbstractVector, β, C::AbstractVector)
    mA, nA = size(A)
    mB = length(B)
    nA == mB || throw(DimensionMismatch("Dimensions must match"))
    length(C) == mA || throw(DimensionMismatch("Dimensions must match"))

    lmul!(β, C)
    (nA == 0 || mB == 0) && return C

    @inbounds for k = 1:mB
        C[k] += α * A.diag[k] * B[k]
    end

    C
end
