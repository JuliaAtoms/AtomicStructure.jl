# * Hartree–Fock operators
# ** Potentials

const HFPotentialOperator{T,B} = RadialOperator{T,B,Diagonal{T,Vector{T}}}

mutable struct HFPotential{kind,O,T,B,OV<:RadialOrbital{T,B},RO<:HFPotentialOperator{T,B},P<:PoissonProblem}
    k::Int
    a::O
    b::O
    av::OV
    bv::OV
    V̂::RO
    poisson::P
end
HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV, V̂::RO, poisson::P) where {O,T,B,OV<:RadialOrbital{T,B},RO<:RadialOperator{T,B},P} =
    HFPotential{kind,O,T,B,OV,RO,P}(k, a, b, av, bv, V̂, poisson)

function HFPotential(kind::Symbol, k::Int, a::O, b::O, av::OV, bv::OV) where {O,T,B,OV<:RadialOrbital{T,B}}
    R = av.mul.factors[1]
    D = Diagonal(Vector{T}(undef, size(R,2)))
    D.diag .= zero(T)
    V̂ = R*D*R'
    poisson = PoissonProblem(k, av, bv, w′=R*D.diag)
    update!(HFPotential(kind, k, a, b, av, bv, V̂, poisson))
end

Base.convert(::Type{HFPotential{kind,O₁,T,B,OV,RO}},
             hfpotential::HFPotential{kind,O₂,T,B,OV,RO,P}) where {kind,O₁,O₂,T,B,OV,RO,P} =
                 HFPotential{kind,O₁,T,B,OV,RO,P}(hfpotential.k,
                                                  hfpotential.a, hfpotential.b, hfpotential.av, hfpotential.bv,
                                                  hfpotential.V̂, hfpotential.poisson)

# *** Direct potential

const DirectPotential{O,T,B,OV,RO} = HFPotential{:direct,O,T,B,OV,RO}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.a), $(Y.b))")

function SCF.update!(p::DirectPotential{O,T,B,OV,RO}; kwargs...) where {O,T,B,OV,RO}
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

mutable struct OrbitalSplitHamiltonian{T, B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       PO<:HFPotentialOperator{T,B},
                                       H<:AbstractOneBodyHamiltonian,
                                       OV<:RadialOrbital{T,B},
                                       Proj}
    R::B
    ĥ::H
    direct_potentials::Vector{Pair{DirectPotential{O,T,B,OV,PO},T}}
    exchange_potentials::Vector{Pair{ExchangePotential{O,T,B,OV,PO},T}}
    projector::Proj
    orbital::O
end

Base.axes(hamiltonian::OrbitalSplitHamiltonian, args...) = axes(hamiltonian.ĥ, args...)
Base.iszero(hamiltonian::OrbitalSplitHamiltonian) =
    iszero(hamiltonian.ĥ) && isempty(hamiltonian.direct_potentials) && isempty(hamiltonian.exchange_potentials)

SCF.update!(hamiltonian::OrbitalSplitHamiltonian; kwargs...) =
    foreach(pc -> update!(pc[1]; kwargs...), hamiltonian.direct_potentials)

function Base.show(io::IO, hamiltonian::OrbitalSplitHamiltonian{T}) where T
    if iszero(hamiltonian)
        write(io, "0")
        return
    end
    multiple_terms = sum([!iszero(hamiltonian.ĥ), !isempty(hamiltonian.direct_potentials),
                          !isempty(hamiltonian.exchange_potentials)]) > 1
    multiple_terms && write(io, "[")
    !iszero(hamiltonian.ĥ) && write(io, "ĥ")
    for (p,c) in hamiltonian.direct_potentials
        s = sign(c)
        write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
    end
    for (p,c) in hamiltonian.exchange_potentials
        s = sign(c)
        write(io, " ", (s < 0 ? "-" : "+"), " $(abs(c))$(p)")
    end
    multiple_terms && write(io, "]")
    write(io, "|")
    show(io, hamiltonian.orbital)
    write(io, "⟩")
end

emptyvec(::V) where {V<:AbstractVector} = V()

function Base.getindex(H::OrbitalSplitHamiltonian, term::Symbol)
    if term == :all
        H
    elseif term == :onebody
        OrbitalSplitHamiltonian(H.R, H.ĥ, emptyvec(H.direct_potentials), emptyvec(H.exchange_potentials), H.projector, H.orbital)
    elseif term == :direct
        OrbitalSplitHamiltonian(H.R, zero(H.ĥ), H.direct_potentials, emptyvec(H.exchange_potentials), H.projector, H.orbital)
    elseif term == :exchange
        OrbitalSplitHamiltonian(H.R, zero(H.ĥ), emptyvec(H.direct_potentials), H.exchange_potentials, H.projector, H.orbital)
    else
        throw(ArgumentError("Unknown Hamiltonian term $term"))
    end
end

# const OrbitalHamiltonian{T,B,O} = Union{OrbitalSplitHamiltonian{T,B,O},RadialOperator{T,B}}

# *** Materialization

const OrbitalSplitHamiltonianMatrixElement{T,V<:AbstractVector,B<:AbstractQuasiMatrix} =
    Mul{<:Tuple,<:Tuple{<:Adjoint{T,V},<:QuasiAdjoint{T,B},<:OrbitalSplitHamiltonian{T,B},B,V}}

const OrbitalSplitHamiltonianMatrixVectorProduct{T,V<:AbstractVector,B<:AbstractQuasiMatrix} =
    Mul{<:Tuple,<:Tuple{<:OrbitalSplitHamiltonian{T,B},B,V}}

Base.eltype(::OrbitalSplitHamiltonianMatrixVectorProduct{T,V,B}) where {T,V,B} = T

function Base.copyto!(dest::RadialOrbital{T,B}, matvec::OrbitalSplitHamiltonianMatrixVectorProduct{T,V,B}) where {T,V,B}
    axes(dest) == axes(matvec) || throw(DimensionMismatch("axes must be the same"))
    R′,v = dest.mul.factors
    A,R,b = matvec.factors

    # One-body operator
    if iszero(A.ĥ)
        v .= zero(T)
    else
        copyto!(v, A.ĥ⋆b)
    end

    for (p,c) in A.direct_potentials
        # # This is how we want to write it, to be basis-agnostic
        # materialize!(MulAdd(c, p.V̂, R*b, one(T), dest))
        materialize!(MulAdd(c, p.V̂.mul.factors[2], b, one(T), dest.mul.factors[2]))
    end
    for (p,c) in A.exchange_potentials
        p.poisson(R*b) # Form exchange potential from conj(p.a)*b
        # # This is how we want to write it, to be basis-agnostic
        # # materialize!(MulAdd(c, p.V̂, p.b, one(T), dest)) # Act on p.b
        materialize!(MulAdd(c, p.V̂.mul.factors[2], p.bv.mul.factors[2], one(T), dest.mul.factors[2])) # Act on p.bv
    end

    projectout!(dest, A.projector)

    dest
end

function Base.similar(matvec::OrbitalSplitHamiltonianMatrixVectorProduct, ::Type{T}) where T
    A,R,b = matvec.factors
    v = similar(b, T)
    R*v
end

LazyArrays.materialize(matvec::OrbitalSplitHamiltonianMatrixVectorProduct{T,V,B}) where {T,V,B} =
    copyto!(similar(matvec, eltype(matvec)), matvec)

function LazyArrays.materialize(matel::OrbitalSplitHamiltonianMatrixElement{T,V,B}) where {T,V,B}
    a,R′,O,R,b = matel.factors
    a*R′*materialize(O⋆R⋆b)
end

# *** Krylov wrapper

function SCF.KrylovWrapper(hamiltonian::OrbitalSplitHamiltonian{T,B,O,PO,LT,OV}) where {T,B,O,PO,LT,OV}
    R = hamiltonian.R
    KrylovWrapper{T,B,OrbitalSplitHamiltonian{T,B,O,PO,LT,OV}}(R, hamiltonian)
end

LinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,B,Hamiltonian}, x::V₂) where {V₁,V₂,T,B,Hamiltonian<:OrbitalSplitHamiltonian} =
    copyto!(A.R*y, A.hamiltonian⋆(A.R*x))
