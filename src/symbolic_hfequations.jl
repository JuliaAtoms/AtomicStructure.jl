# This manner of constructing the Hartree–Fock equations, using
# symbolic methods, is illustrative, but slow.

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

    ĥ = one_body_hamiltonian(atom, orbital)

    eq_terms = equation_terms(equation)
    korb = Ket(orbital)

    OV = typeof(view(atom, 1))
    PO = typeof(R*Diagonal(Vector{T}(undef, size(R,2)))*R')
    direct_potentials = Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}()
    exchange_potentials = Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}()

    for t ∈ eq_terms
        args = t.args
        if isproportional(t, LagrangeMultiplier)
            lm = findfactor(t, LagrangeMultiplier)
            # TODO: Apparently, off-diagonal Lagrange multipliers are
            # only necessary between orbitals of the same symmetry,
            # but which are varied separately (or only one of them is
            # varied, e.g. when 1s is frozen and 2s is varied), and
            # the other orbital needs to be projected out at each
            # step. /However/ the rotation may be necessary in any
            # case, to speed up convergence, since there are
            # infinitely many solutions.
            isdiagonal(lm) ||
                throw(ArgumentError("Orbital rotation necessary of off-diagonal Lagrange multipliers ($(lm)) not yet implemented."))
        elseif isproportional(t, Sym(:𝓛))
            isproportional(t, korb) ||
                throw(ArgumentError("One-body term $(t) not pertaining to orbital under consideration ($(orbital))"))
            cL̂ = -convert(T,(t/korb)/Sym(:𝓛))
            @assert cL̂ ≈ one(T)
        elseif isproportional(t, SlaterPotential)
            Y = findfactor(t, SlaterPotential)
            other = Y.a == orbital ? Y.b : Y.a
            kind = isexchange(Y, orbital) && other != orbital ? :exchange : :direct
            k = log(Y)
            v = view(atom, other)
            # We negate the coefficients of the Slater potentials here
            # such that ĥ has +1 as coefficient. Furthermore, we
            # divide by 2, since we want the HF equation to pertain to
            # singly occupied orbital, i.e. we use ĥ as the one-body
            # operator, not 𝓛 = 2ĥ. This has no effect on HF
            # iterations, since we want to converge [H-E]|orb⟩ = 0 –
            # which is not affected by this choice (E then pertains to
            # the removal of a single electron from the orbital) – but
            # it is useful for time propagation later.
            c = -convert(T, (t*Sym(:r)/Y)/(kind == :direct ? korb : Ket(other)))/2

            pc = HFPotential(kind, k, other, v) => c
            if kind == :direct
                push!(direct_potentials, pc)
            else
                push!(exchange_potentials, pc)
            end
        else
            throw(ArgumentError("Unknown Hartree–Fock equation term $(t)"))
        end
    end

    hamiltonian = OrbitalSplitHamiltonian(ĥ, direct_potentials, exchange_potentials)

    HFEquation(atom, equation, orbital, view(atom, orbital), hamiltonian)
end
