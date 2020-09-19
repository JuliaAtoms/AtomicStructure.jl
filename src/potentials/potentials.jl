abstract type AbstractPotential{T} end

include("table_of_elements.jl")

# * Point charge nucleus

struct PointCharge{T} <: AbstractPotential{T}
    Z::T
end

PointCharge(s::Symbol) = PointCharge(element_number(s))

macro pc_str(s)
    :(PointCharge(Symbol($s)))
end

function Base.show(io::IO, p::PointCharge{I}) where {I<:Integer}
    element = table_of_elements[p.Z]
    write(io, "Z = $(p.Z) [$(element[2][1])]")
end

Base.show(io::IO, p::PointCharge) =
    write(io, "Z = $(p.Z)")

(p::PointCharge{T})( ::O, r::U) where {T,O,U} = -p.Z/r
(p::PointCharge{T})(o::O, r::VU) where {T,O,U,VU<:AbstractVector{U}} = p.(Ref(o), r)

charge(p::PointCharge) = p.Z

ground_state(p::PointCharge{<:Integer}) =
    table_of_elements[p.Z][2][2]

ground_state(p::PointCharge) =
    throw(ArgumentError("Ground state configuration for nuclear charge of Z = $(p.Z) unknown"))

islocal(::PointCharge) = true

# * Yukawa

struct Yukawa{T} <: AbstractPotential{T}
    g::T
    α::T
    m::T
end

(p::Yukawa)( ::O, r::U) where {O,U} = -p.g^2*exp(-p.α*p.m*r)/r
(p::Yukawa)(o::O, r::VU) where {O,VU<:AbstractVector} = p.(Ref(o), r)

islocal(::Yukawa) = true

# * ECPs

abstract type AbstractEffectiveCorePotential{T} <: AbstractPotential{T} end

islocal(::AbstractEffectiveCorePotential) = false

include("gaussian_expansions.jl")

include("nonrelativistic_effective_core_potentials.jl")
include("relativistic_effective_core_potentials.jl")
const FullOrScalarRelativisticEffectiveCorePotential = Union{RelativisticEffectiveCorePotential,ScalarSORelativisticEffectiveCorePotential}

include("parse.jl")
include("misc_effective_core_potentials.jl")

# * Exports

export AbstractPotential,
    PointCharge, Yukawa,
    EffectiveCorePotential, RelativisticEffectiveCorePotential,
    @pc_str, @ECP_str,
    element_number,
    charge, ground_state, islocal,
    spin_orbit_potential
