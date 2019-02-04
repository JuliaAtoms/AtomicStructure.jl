# * Generation of one-body Hamiltonian

function one_body_hamiltonian(::Type{Tuple}, atom::Atom{T,T,B,O₁,TC,C,CM}, orb::O₂) where {T,B,O₁,TC,C,CM,O₂}
    R = radial_basis(atom)

    D = Derivative(axes(R,1))
    Tm = R'D'D*R
    Tm /= -2

    ℓ = getℓ(orb)

    V = Matrix(r -> ℓ*(ℓ+1)/(2r^2) + atom.potential(orb, r), R)

    R*Tm*R', R*V*R'
end
one_body_hamiltonian(atom, orb) = +(one_body_hamiltonian(Tuple, atom, orb)...)

# * OneBodyHamiltonian types

abstract type AbstractOneBodyHamiltonian end
Base.iszero(::AbstractOneBodyHamiltonian) = false

# ** ZeroOneBodyHamiltonian

struct ZeroOneBodyHamiltonian{M,N} <: AbstractOneBodyHamiltonian
    a::Tuple{M,N}
end

Base.axes(ĥ::ZeroOneBodyHamiltonian) = ĥ.a
Base.axes(ĥ::ZeroOneBodyHamiltonian,i) = ĥ.a[i]

Base.iszero(::ZeroOneBodyHamiltonian) = true

Base.show(io::IO, ĥ::ZeroOneBodyHamiltonian) =
    write(io, "0")

# ** AtomicOneBodyHamiltonian
struct AtomicOneBodyHamiltonian{LT,O} <: AbstractOneBodyHamiltonian
    op::LT
    orbital::O
end

AtomicOneBodyHamiltonian(atom::Atom, orb::AbstractOrbital) =
    AtomicOneBodyHamiltonian(one_body_hamiltonian(atom, orb), orb)

Base.axes(ĥ::AtomicOneBodyHamiltonian, args...) =
    axes(ĥ.op, args...)

Base.zero(ĥ::AtomicOneBodyHamiltonian) =
    ZeroOneBodyHamiltonian(axes(ĥ))

LazyArrays.:(⋆)(ĥ::AtomicOneBodyHamiltonian, ϕ::AbstractVector) =
    ĥ.op.mul.factors[2] ⋆ ϕ

Base.show(io::IO, ĥ::AtomicOneBodyHamiltonian) =
    write(io, "ĥ($(ĥ.orbital))")

# * Diagonalization of one-body Hamiltonian

struct ShiftInvert{M}
    A⁻¹::M
end

Base.size(S::ShiftInvert, args...) = size(S.A⁻¹, args...)
Base.eltype(S::ShiftInvert) = eltype(S.A⁻¹)

LinearAlgebra.mul!(y, S::ShiftInvert, x) =
    ldiv!(y, S.A⁻¹, x)

function diagonalize_one_body(H::RadialOperator, nev::Int;
                              method::Symbol=:arnoldi_shift_invert, tol=1e-10, σ=-1,
                              verbosity=0, io=stdout)
    Hm = matrix(H)
    verbosity > 2 && println(io, "Diagonalizing via $(method)")
    if method == :arnoldi || method == :arnoldi_shift_invert
        A,target = if method == :arnoldi
            Hm,SR()
        else
            ShiftInvert(factorize(Hm - σ*I)),LR()
        end
        schur,history = partialschur(A, nev=nev, tol=tol, which=target)
        verbosity > 3 && println(io, history)
        λ = if method == :arnoldi
            diag(schur.R)
        else
            θ = diag(schur.R)
            verbosity > 2 && println(io, "Schur values: $θ")
            σ .+ inv.(θ)
        end
        length(λ) < nev &&
            error("Could not converge the requested orbitals: $(history)")
        λ,schur.Q
    elseif method == :eigen
        ee = eigen(Hm)
        ee.values[1:nev],ee.vectors[:,1:nev]
    else
        throw(ArgumentError("Unknown diagonalization method $(method)"))
    end
end
