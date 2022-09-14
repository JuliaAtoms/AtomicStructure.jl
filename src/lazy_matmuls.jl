for MT in [:Diagonal, :SymTridiagonal, :Tridiagonal]
    @eval begin
        function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, <:Any,
                                                    <:RadialOperator{<:Any, <:$MT},
                                                    <:RadialOrbital, <:RadialOrbital})
            mul!(ma.C.args[2], ma.A.args[2], ma.B.args[2], ma.α, ma.β)
            ma.C
        end

        function LazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, <:Any,
                                                    <:RadialOperator{<:Any, Adjoint{<:Any, <:$MT}},
                                                    <:RadialOrbital, <:RadialOrbital})
            mul!(ma.C.args[2], ma.A.args[2], ma.B.args[2], ma.α, ma.β)
            ma.C
        end
    end
end
