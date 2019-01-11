# * Hartree–Fock operators

mutable struct HFPotential{kind,O,T,B,V<:RadialOrbital{T,B}}
    k::Int
    orbital::O
    v::V
end
HFPotential(kind::Symbol, k::Int, orbital::O, v::V) where {O,T,B,V<:RadialOrbital{T,B}} =
    HFPotential{kind,O,T,B,V}(k, orbital, v)

const DirectPotential{O,T,B,V} = HFPotential{:direct,O,T,B,V}

Base.show(io::IO, Y::DirectPotential) =
    write(io, "r⁻¹×Y", to_superscript(Y.k), "($(Y.orbital), $(Y.orbital))")

const ExchangePotential{O,T,B,V} = HFPotential{:exchange,O,T,B,V}

Base.show(io::IO, Y::ExchangePotential) =
    write(io, "|$(Y.orbital)⟩r⁻¹×Y", to_superscript(Y.k), "($(Y.orbital), ●)")

mutable struct OrbitalSplitHamiltonian{T,ΦT, #<:RadialCoeff{T},
                                       B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       # RO<:RadialOperator{ΦT,B},
                                       LT, #<:Union{RO,NTuple{<:Any,RO}},
                                       V<:RadialOrbital{ΦT,B}}
    L̂::LT
    cL̂::T
    direct_potentials::Vector{Pair{DirectPotential{O,ΦT,B,V},T}}
    exchange_potentials::Vector{Pair{ExchangePotential{O,ΦT,B,V},T}}
end

function Base.show(io::IO, hamiltonian::OrbitalSplitHamiltonian{T}) where T
    hamiltonian.cL̂ != one(T) && show(io, hamiltonian.cL̂)
    write(io, "𝓛")
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

mutable struct HFEquation{T, ΦT, #<:RadialCoeff{T},
                          B<:AbstractQuasiMatrix,
                          O<:AbstractOrbital,
                          A<:Atom{T,ΦT,B,O},
                          E<:Symbolic,
                          M# <:OrbitalHamiltonian{T,ΦT,B,O}
                          }
    atom::A
    equation::E
    orbital::O
    ϕ::RadialOrbital{T,B}
    hamiltonian::M
end

equation_terms(equation::Number) = [equation]
equation_terms(equation::SymExpr) = equation.op == Sym(:+) ? equation.args : [equation.args]

isproportional(equation::Number, ::Type{T}) where T = false
isproportional(equation::SymExpr, ::Type{T}) where T =
    equation.op == Sym(:*) && any(a isa T for a in equation.args)

isproportional(equation::Number, s::Number) =
    equation == s
isproportional(equation::SymExpr, s::Number) =
    equation.op == Sym(:*) && s ∈ equation.args

findfactor(equation::SymExpr, ::Type{T}) where T =
    equation.args[findfirst(a -> a isa T, equation.args)]

function HFEquation(atom::A, equation::E, orbital::O) where {T,ΦT, #<:RadialCoeff{T},
                                                             B<:AbstractQuasiMatrix,
                                                             O<:AbstractOrbital,
                                                             A<:Atom{T,ΦT,B,O},
                                                             E<:Symbolic}
    R = radial_basis(atom)

    L̂ = one_body_hamiltonian(atom, orbital)
    cL̂ = one(T)

    eq_terms = equation_terms(equation)
    korb = Ket(orbital)

    direct_potentials,exchange_potentials = map([:direct,:exchange]) do kind
        pred = t -> (isproportional(t, SlaterPotential) &&
                     ((kind == :direct) ⊻ isexchange(findfactor(t, SlaterPotential), orbital)))
        map(filter(pred, eq_terms)) do t
            Y = findfactor(t, SlaterPotential)
            k = log(Y)
            other = Y.a == orbital ? Y.b : Y.a
            v = view(atom, other)
            c = convert(T, (t*Sym(:r)/Y)/(kind == :direct ? korb : Ket(other)))

            HFPotential(kind, k, other, v) => c
        end
    end

    for t ∈ eq_terms
        args = t.args
        if isproportional(t, LagrangeMultiplier)
            lm = findfactor(t, LagrangeMultiplier)
            isdiagonal(lm) ||
                throw(ArgumentError("Orbital rotation necessary of off-diagonal Lagrange multipliers ($(lm)) not yet implemented."))
        elseif isproportional(t, Sym(:𝓛))
            isproportional(t, korb) ||
                throw(ArgumentError("One-body term $(t) not pertaining to orbital under consideration ($(orbital))"))
            cL̂ = convert(T,(t/korb)/Sym(:𝓛))
        elseif isproportional(t, SlaterPotential)

        else
            throw(ArgumentError("Unknown Hartree–Fock equation term $(t)"))
        end
    end

    hamiltonian = OrbitalSplitHamiltonian(L̂, cL̂, direct_potentials, exchange_potentials)

    HFEquation(atom, equation, orbital, view(atom, orbital), hamiltonian)
end

# # We only define the action of the Hamiltonian on an arbitrary
# # vector. The solution of the equation is implemented via iterative
# # eigenfactorization.
# LinearAlgebra.mul!(y, hfeq::HFEquation, x) = mul!(y, hfeq.hamiltonian, x)
# function Base.:(*)(hfeq::HFEquation, x)
#     # y = similar(x)
#     # y .= hfeq.hamiltonian ⋆ x
#     # y
# end

energy(hfeq::HFEquation{E,O,M}) where {E,O,M} =
    Inf # materialize(hfeq.ϕ' ⋆ hfeq.hamiltonian ⋆ hfeq.ϕ)[1]

function Base.show(io::IO, hfeq::HFEquation)
    write(io, "Hartree–Fock equation: 0 = [𝓗  - E($(hfeq.orbital))]|$(hfeq.orbital)⟩ = [$(hfeq.hamiltonian) - E($(hfeq.orbital))]|$(hfeq.orbital)⟩")
    write(io, "\n    ⟨$(hfeq.orbital)| 𝓗 |$(hfeq.orbital)⟩ = $(energy(hfeq))")
end

# * Setup Hartree–Fock equations

"""
    hf_equations(csf[; verbosity=0])

Derive the Hartree–Fock equations for the non-relativistic
configuration state function `csf`. All constituent orbitals are
assumed spectroscopic, i.e. they are all assigned Lagrange multipliers
to ensure their orthonormality. Additionally, orbitals of the same `ℓ`
but different `n` are assigned off-diagonal Lagrange multipliers.
"""
function hf_equations(atom::Atom, csf::NonRelativisticCSF; verbosity=0)
    verbosity > 0 &&
        println("Deriving HF equations for $(csf)")

    orbitals = csf.config.orbitals
    # energy_expression has to be provided by an angular momentum
    # library.
    eng = energy_expression(csf; verbosity=verbosity-3)[2][1]
    λs = lagrange_multipliers(orbitals)
    eng += λs

    if verbosity > 2
        println("Orbitals: $(orbitals)")
        println("Lagrange multipliers: $(λs)")
    end
    verbosity > 1 && println("E = $eng\n")

    map(peel(csf.config)) do (orb,occ,state)
        hf_eq = diff(eng, orb, occ)
        verbosity > 2 && println("∂[$(orb)] E = $(hf_eq)")
        HFEquation(atom, hf_eq, orb)
    end
end

function Base.diff(atom::Atom; kwargs...)
    length(atom.csfs) > 1 &&
        throw(ArgumentError("Cannot derive Hartree–Fock equations for a multi-configurational atom"))

    hf_equations(atom, atom.csfs[1]; kwargs...)
end

Base.view(fock::Fock{A,E}, args...) where {A<:Atom,E} =
    view(fock.quantum_system.radial_orbitals.mul.factors[2], args...)

# * Solve one iteration of the Hartree–Fock equations

function LinearAlgebra.ldiv!(fock::Fock{A,E}, c::M) where {A<:Atom,E,M<:AbstractMatrix}
    c .= rand() # Bogus solution of secular equations
end
