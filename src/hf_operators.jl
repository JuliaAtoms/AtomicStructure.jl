# * Hartree–Fock operators
# ** Potentials

const HFPotentialOperator{T,B} = RadialOperator{T,B,Diagonal{T,Vector{T}}}

mutable struct HFPotential{kind,O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B},RO<:HFPotentialOperator{T,B},P<:PoissonProblem}
    k::Int
    a::O
    b::O
    av::OV
    bv::OV
    V̂::RO
    poisson::P
end
HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV, V̂::RO, poisson::P) where {O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B},RO<:RadialOperator{T,B},P} =
    HFPotential{kind,O,T,ΦT,B,OV,RO,P}(k, a, b, av, bv, V̂, poisson)

function HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV) where {O,T,ΦT<:RadialCoeff{T},B,OV<:RadialOrbital{ΦT,B}}
    R = av.mul.factors[1]
    D = Diagonal(Vector{T}(undef, size(R,2)))
    D.diag .= zero(T)
    V̂ = R*D*R'
    poisson = PoissonProblem(k, av, bv, w′=R*D.diag)
    update!(HFPotential(kind, k, a, b, av, bv, V̂, poisson))
end

Base.convert(::Type{HFPotential{kind,O₁,T,ΦT,B,OV,RO}},
             hfpotential::HFPotential{kind,O₂,T,ΦT,B,OV,RO,P}) where {kind,O₁,O₂,T,ΦT,B,OV,RO,P} =
                 HFPotential{kind,O₁,T,ΦT,B,OV,RO,P}(hfpotential.k,
                                                     hfpotential.a, hfpotential.b, hfpotential.av, hfpotential.bv,
                                                     hfpotential.V̂, hfpotential.poisson)

# *** Direct potential

const DirectPotential{O,T,ΦT,B,OV,RO} = HFPotential{:direct,O,T,ΦT,B,OV,RO}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.a), $(Y.b))")

function SCF.update!(p::DirectPotential{O,T,ΦT,B,OV,RO}; kwargs...) where {O,T,ΦT,B,OV,RO}
    p.poisson(;kwargs...)
    p
end

# *** Exchange potential

const ExchangePotential{O,T,B,OV,RO} = HFPotential{:exchange,O,T,B,OV,RO}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.b)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.a), ●)")

# We can't update the exchange potentials, since they depend on the
# orbital they act on.
SCF.update!(p::ExchangePotential; kwargs...) = p

# ** Split Hamiltonian

mutable struct OrbitalSplitHamiltonian{T,ΦT, #<:RadialCoeff{T},
                                       B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       PO<:HFPotentialOperator{T,B},
                                       LT, #<:Union{RO,NTuple{<:Any,RO}},
                                       OV<:RadialOrbital{ΦT,B},
                                       Proj}
    ĥ::LT
    direct_potentials::Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}
    exchange_potentials::Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}
    projector::Proj
end

Base.axes(hamiltonian::OrbitalSplitHamiltonian, args...) = axes(hamiltonian.ĥ, args...)

SCF.update!(hamiltonian::OrbitalSplitHamiltonian; kwargs...) =
    foreach(pc -> update!(pc[1]; kwargs...), hamiltonian.direct_potentials)

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

emptyvec(::V) where {V<:AbstractVector} = V()

function Base.getindex(H::OrbitalSplitHamiltonian, term::Symbol)
    if term == :all
        H
    elseif term == :onebody
        OrbitalSplitHamiltonian(H.ĥ, emptyvec(H.direct_potentials), emptyvec(H.exchange_potentials), H.projector)
    elseif term == :direct
        OrbitalSplitHamiltonian(zero(H.ĥ), H.direct_potentials, emptyvec(H.exchange_potentials), H.projector)
    elseif term == :exchange
        OrbitalSplitHamiltonian(zero(H.ĥ), emptyvec(H.direct_potentials), H.exchange_potentials, H.projector)
    else
        throw(ArgumentError("Unknown Hamiltonian term $term"))
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
        # # This is how we want to write it, to be basis-agnostic
        # materialize!(MulAdd(c, p.V̂, R*b, one(T), dest))
        materialize!(MulAdd(c, p.V̂.mul.factors[2], b, one(T), dest.mul.factors[2]))
    end
    for (p,c) in A.exchange_potentials
        p.poisson(R*b) # Form exchange potential from conj(p.a)*b
        # # This is how we want to write it, to be basis-agnostic
        # # materialize!(MulAdd(c, p.V̂, p.b, one(T), dest)) # Act on p.b
        materialize!(MulAdd(c, p.V̂.mul.factors[2], p.b.mul.factors[2], one(T), dest.mul.factors[2])) # Act on p.b
    end

    projectout!(dest, A.projector)

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

function SCF.KrylovWrapper(hamiltonian::OrbitalSplitHamiltonian{T,ΦT,B,O,PO,LT,OV}) where {T,ΦT,B,O,PO,LT,OV}
    R = hamiltonian.ĥ.mul.factors[1]
    KrylovWrapper{ΦT,B,OrbitalSplitHamiltonian{T,ΦT,B,O,PO,LT,OV}}(R, hamiltonian)
end

LinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,B,Hamiltonian}, x::V₂) where {V₁,V₂,T,B,Hamiltonian<:OrbitalSplitHamiltonian} =
    copyto!(A.R*y, A.hamiltonian⋆(A.R*x))
