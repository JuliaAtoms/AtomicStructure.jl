# * Hartree–Fock equation

mutable struct HFEquation{T, ΦT, #<:RadialCoeff{T},
                          B<:AbstractQuasiMatrix,
                          O<:AbstractOrbital,
                          A<:Atom{T,ΦT,B,O},
                          E,
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

function HFEquation(atom::A, (one_body,(two_body,multipole_terms))::E,
                    orbital::O) where {T,ΦT, #<:RadialCoeff{T},
                                       B<:AbstractQuasiMatrix,
                                       O<:AbstractOrbital,
                                       A<:Atom{T,ΦT,B,O},
                                       E<:Tuple{Vector{<:OneBodyHamiltonian},Tuple{Vector{<:DirectExchangePotentials},Vector{NTuple{2,Vector{Pair{Int,Float64}}}}}}}
    R = radial_basis(atom)

    ĥ = one_body_hamiltonian(atom, orbital)

    korb = Ket(orbital)

    OV = typeof(view(atom, 1))
    PO = typeof(R*Diagonal(Vector{T}(undef, size(R,2)))*R')
    direct_potentials = Vector{Pair{DirectPotential{O,T,ΦT,B,OV,PO},T}}()
    exchange_potentials = Vector{Pair{ExchangePotential{O,T,ΦT,B,OV,PO},T}}()

    count(iszero.(one_body)) == 1 ||
        throw(ArgumentError("There can only be one one-body Hamiltonian per orbital"))

    all(AngularMomentumAlgebra.isdiagonal.(two_body)) ||
        throw(ArgumentError("Non-diagonal repulsion potentials not yet supported"))

    for (tb,mpt) ∈ zip(two_body,multipole_terms)
        tb.o.v == orbital ||
            throw(ArgumentError("Repulsion potential $(tb) not pertaining to $(orbital)"))

        other = tb.a
        v = view(atom, other)

        for (k,c) in mpt[1]
            push!(direct_potentials,
                  HFPotential(:direct, k, other, v) => c)
        end
        for (k,c) in mpt[2]
            push!(exchange_potentials,
                  HFPotential(:exchange, k, other, v) => c)
        end
    end

    hamiltonian = OrbitalSplitHamiltonian(ĥ, direct_potentials, exchange_potentials)

    HFEquation(atom, (one_body,two_body), orbital, view(atom, orbital), hamiltonian)
end

SCF.energy(hfeq::HFEquation{E,O,M}) where {E,O,M} = (hfeq.ϕ' * hfeq.hamiltonian * hfeq.ϕ)[1]

function Base.show(io::IO, hfeq::HFEquation)
    write(io, "Hartree–Fock equation: 0 = [𝓗  - E($(hfeq.orbital))]|$(hfeq.orbital)⟩ = [$(hfeq.hamiltonian) - E($(hfeq.orbital))]|$(hfeq.orbital)⟩")
    EHa = SCF.energy(hfeq)
    write(io, "\n    ⟨$(hfeq.orbital)| 𝓗 |$(hfeq.orbital)⟩ = $(EHa) Ha = $(27.211EHa) eV")
end

# * Setup Hartree–Fock equations

function hf_equations(csf::NonRelativisticCSF, eng::Number; verbosity=0)
    pconfig = peel(csf.config)
    orbitals = pconfig.orbitals

    λs = lagrange_multipliers(orbitals)
    eng += λs

    if verbosity > 2
        println("Orbitals: $(orbitals)")
        println("Lagrange multipliers: $(λs)")
    end
    verbosity > 1 && println("E = $eng\n")

    map(pconfig) do (orb,occ,state)
        equation = diff(eng, orb, occ)
        verbosity > 2 && println("\n∂[$(orb)] E = $(equation)")
        orb => equation
    end
end

"""
    hf_equations(csf[; verbosity=0])

Derive the Hartree–Fock equations for the non-relativistic
configuration state function `csf`. All constituent orbitals are
assumed spectroscopic, i.e. they are all assigned Lagrange multipliers
to ensure their orthonormality. Additionally, orbitals of the same `ℓ`
but different `n` are assigned off-diagonal Lagrange multipliers.
"""
function hf_equations(csf::NonRelativisticCSF; verbosity=0)
    verbosity > 0 &&
        println("Deriving HF equations for $(csf)")

    # energy_expression has to be provided by an angular momentum
    # library.
    eng = energy_expression(csf; verbosity=verbosity-3)[2][1]
    hf_equations(csf, eng; verbosity=verbosity)
end

function hf_equations(config::Configuration{O}; verbosity=0) where {O<:SpinOrbital}
    verbosity > 0 &&
        println("Deriving HF equations for $(config)")
    h = one_body_hamiltonian_matrix(O, [config])[1]
    HC = two_body_hamiltonian_matrix(O, [config])[1]
    map(config.orbitals) do orb
        corb = Conjugate(orb)
        multipole_terms = multipole_expand.(HC.integrals)
        orb,(diff.(h.integrals, Ref(corb)), (diff.(HC.integrals, Ref(corb)), multipole_terms))
    end
end

"""
    diff(atom)

Differentiate the energy expression associated with the `atom`'s
CSF(s) with respect to the atomic orbitals to derive the Hartree–Fock
equations for the orbitals.
"""
function Base.diff(atom::Atom; verbosity=0)
    length(atom.configurations) > 1 &&
        throw(ArgumentError("Cannot derive Hartree–Fock equations for a multi-configurational atom"))

    map(hf_equations(atom.configurations[1]; verbosity=verbosity)) do (orb,equation)
        verbosity > 3 && display(equation)
        hfeq = HFEquation(atom, equation, orb)
        verbosity > 2 && println(hfeq)
        hfeq
    end
end

Base.view(fock::Fock{A,E}, args...) where {A<:Atom,E} =
    view(fock.quantum_system.radial_orbitals.mul.factors[2], args...)

# * Solve one iteration of the Hartree–Fock equations

function LinearAlgebra.ldiv!(fock::Fock{A,E}, c::M;
                             verbosity=0, method=:arnoldi,
                             tol=1e-10, kwargs...) where {A<:Atom,E,M<:AbstractVecOrMat}
    verbosity > 0 && println("Solving secular equations using $(method)")
    atom = fock.quantum_system
    R = radial_basis(atom)
    for (j,eq) in enumerate(fock.equations)
        print_block() do io
            # Ideally, all direct potentials should be shared among the
            # equations, such that they are only calculated once per
            # iteration.
            verbosity > 1 && println(io, eq)
            update!(eq.hamiltonian; verbosity=max(0,verbosity-2), io=io)
            if method==:arnoldi
                K = KrylovWrapper(eq.hamiltonian)
                schur,history = partialschur(K, nev=1, tol=tol, which=SR())
                vcj = view(c,:,j)
                copyto!(vcj, schur.Q[:,1])
                norm_rot!(R*vcj)

                verbosity > 2 && println(io, "Secular equation solution: ", history)
                verbosity > 2 && println(io, "Change in equation $j: ",
                                        norm(c[:,j] - eq.ϕ.mul.factors[2]))
            else
                throw(ArgumentError("Unknown diagonalization method $(method)"))
            end
        end
    end
    fock
end
