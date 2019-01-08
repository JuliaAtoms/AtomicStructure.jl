abstract type AbstractPotential{T} end

# * Point charge nucleus

struct PointCharge{T} <: AbstractPotential{T}
    Z::T
end

PointCharge(s::Symbol) = PointCharge(element_number(s))

macro pc_str(s)
    :(PointCharge(Symbol($s)))
end

Base.show(io::IO, p::PointCharge{I}) where {I<:Integer} =
    write(io, "Z = $(p.Z) [$(table_of_elements[p.Z])]")

Base.show(io::IO, p::PointCharge) =
    write(io, "Z = $(p.Z)")

(p::PointCharge{T})(ℓ::Int, r::U) where {T,U} = -p.Z/r

charge(p::PointCharge) = p.Z

export PointCharge, @pc_str, charge
