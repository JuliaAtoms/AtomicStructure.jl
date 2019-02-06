struct Projector{T,B<:AbstractQuasiMatrix,RO<:RadialOrbital{T,B}}
    ϕs::Vector{RO}
end

projectout!(y::RO, ::Nothing) where RO = y

function projectout!(y::RO, projector::Proj) where {RO,Proj<:Projector}
    yc = y.mul.factors[2]
    
    for ϕ in projector.ϕs
        c = ϕ'y
        yc .-= c*ϕ.mul.factors[2]
        # y -= c*ϕ
    end
    y
end
