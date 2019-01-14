# * Hartree–Fock operators
# ** Potentials

const HFPotentialOperator{T,B} = RadialOperator{T,B,Diagonal{T,Vector{T}}}

mutable struct HFPotential{kind,O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B},RO<:HFPotentialOperator{T,B}}
    k::Int
    orbital::O
    ov::OV
    V̂::RO
end
HFPotential(kind::Symbol, k::Int, orbital::O, ov::OV, V̂::RO) where {O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B},RO<:RadialOperator{T,B}} =
    HFPotential{kind,O,T,ΦT,B,OV,RO}(k, orbital, ov, V̂)
    
function HFPotential(kind::Symbol, k::Int, orbital::O, ov::OV) where {O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B}}
    R = ov.mul.factors[1]
    V̂ = R*Diagonal(Vector{T}(undef, size(R,2)))*R'
    update!(HFPotential(kind, k, orbital, ov, V̂))
end

# *** Direct potential

const DirectPotential{O,T,ΦT,B,OV,RO} = HFPotential{:direct,O,T,ΦT,B,OV,RO}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.orbital), $(Y.orbital))")

function update!(p::DirectPotential{O,T,ΦT,B,OV,RO}) where {O,T,ΦT,B,OV,RO}
    # TODO: Solve Poisson
    p.V̂.mul.factors[2].diag .= zero(T)
    p
end

# *** Exchange potential

const ExchangePotential{O,T,B,OV,RO} = HFPotential{:exchange,O,T,B,OV,RO}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.orbital)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.orbital), ●)")

update!(p::ExchangePotential) = p

# ** Krylov wrapper
mutable struct KrylovWrapper{T,B<:AbstractQuasiMatrix,Hamiltonian}
    R::B
    hamiltonian::Hamiltonian
end

Base.eltype(A::KrylovWrapper{T}) where T = T
Base.size(A::KrylovWrapper) = (size(A.R,2),size(A.R,2))
Base.size(A::KrylovWrapper, i) = size(A)[i]

function Base.show(io::IO, kw::KrylovWrapper{T,B,Hamiltonian}) where {T,B,Hamiltonian}
    write(io, "KrylovWrapper{$T} of size $(size(kw))")
end

LinearAlgebra.mul!(y::V, A::KrylovWrapper{T,B,Hamiltonian}, x::V) where {V,T,B,Hamiltonian} =
    copyto!(A.R*y, A.hamiltonian⋆(A.R*x))

# ** Split Hamiltonian

mutable struct OrbitalSplitHamiltonian{T,ΦT, #<:RadialCoeff{T},
                                       B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       PO<:HFPotentialOperator{T,B},
                                       LT, #<:Union{RO,NTuple{<:Any,RO}},
                                       OV<:RadialOrbital{ΦT,B}}
    ĥ::LT
    direct_potentials::Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}
    exchange_potentials::Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}
end

Base.axes(hamiltonian::OrbitalSplitHamiltonian, args...) = axes(hamiltonian.ĥ, args...)

update!(hamiltonian::OrbitalSplitHamiltonian) = 
    foreach(pc -> update!(pc[1]), hamiltonian.direct_potentials)

function Base.show(io::IO, hamiltonian::OrbitalSplitHamiltonian{T}) where T
    write(io, "ĥ")
    for (p,c) in hamiltonian.direct_potentials
        s = sign(c)
        write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
    end
    for (p,c) in hamiltonian.exchange_potentials
        s = sign(c)
        write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
    end
end

# const OrbitalHamiltonian{T,ΦT,B,O} = Union{OrbitalSplitHamiltonian{T,ΦT,B,O},RadialOperator{T,B}}

# *** Materialization

const OrbitalSplitHamiltonianMatrixElement{T,ΦT<:RadialCoeff{T},V<:AbstractVector,B<:AbstractQuasiMatrix} =
    Mul{<:Tuple,<:Tuple{<:Adjoint{T,V},<:QuasiAdjoint{T,B},<:OrbitalSplitHamiltonian{T,ΦT,B},B,V}}

const OrbitalSplitHamiltonianMatrixVectorProduct{T,ΦT<:RadialCoeff{T},V<:AbstractVector,B<:AbstractQuasiMatrix} =
    Mul{<:Tuple,<:Tuple{<:OrbitalSplitHamiltonian{T,ΦT,B},B,V}}

Base.eltype(::OrbitalSplitHamiltonianMatrixVectorProduct{T,ΦT,V,B}) where {T,ΦT,V,B} = ΦT

function Base.copyto!(dest::RadialOrbital{ΦT,B}, matvec::OrbitalSplitHamiltonianMatrixVectorProduct{T,ΦT,V,B}) where {T,ΦT,V,B}
    axes(dest) == axes(matvec) || throw(DimensionMismatch("axes must be the same"))
    R′,v = dest.mul.factors
    A,R,b = matvec.factors

    # One-body operator
    ĥ = A.ĥ.mul.factors[2]
    copyto!(v, ĥ⋆b)

    for (p,c) in A.direct_potentials
    end
    for (p,c) in A.exchange_potentials
    end
    
    dest
end

function Base.similar(matvec::OrbitalSplitHamiltonianMatrixVectorProduct, ::Type{T}) where T
    A,R,b = matvec.factors
    v = similar(b, T)
    R*v
end

LazyArrays.materialize(matvec::OrbitalSplitHamiltonianMatrixVectorProduct{T,ΦT,V,B}) where {T,ΦT,V,B} = copyto!(similar(matvec, eltype(matvec)), matvec)

function LazyArrays.materialize(matel::OrbitalSplitHamiltonianMatrixElement{T,ΦT,V,B}) where {T,ΦT,V,B}
    a,R′,O,R,b = matel.factors
    a*R′*materialize(O⋆R⋆b)
end

# *** Krylov wrapper

function KrylovWrapper(hamiltonian::OrbitalSplitHamiltonian{T,ΦT,B,O,PO,LT,OV}) where {T,ΦT,B,O,PO,LT,OV}
    R = hamiltonian.ĥ.mul.factors[1]
    KrylovWrapper{ΦT,B,OrbitalSplitHamiltonian{T,ΦT,B,O,PO,LT,OV}}(R, hamiltonian)
end
