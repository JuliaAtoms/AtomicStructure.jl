var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Atoms.jl-1",
    "page": "Home",
    "title": "Atoms.jl",
    "category": "section",
    "text": "Documentation for Atoms.jl"
},

{
    "location": "density_matrices/#",
    "page": "Density matrices",
    "title": "Density matrices",
    "category": "page",
    "text": ""
},

{
    "location": "density_matrices/#Density-matrices-1",
    "page": "Density matrices",
    "title": "Density matrices",
    "category": "section",
    "text": "The density matrix operator is defined asbeginequation\nhatrho_ij defd ketbraij\nendequationIf the wavefunction is approximated using Slater determinants, where one Slater determinant of N electrons is defined asbeginequation\nPhi(123N) =\nfrac1sqrtN\nleftbeginmatrix\nchi_1(1)\nchi_2(2)\n\nchi_N(N)\nendmatrixright\nendequationwhere 123N denote the N different electronic coordinates (spatial and spin) and chi_i are the N different spin-orbitals, the transition density matrix for all N electrons, between two Slater determinants Phi_A and Phi_B is given bybeginequation\nbeginaligned\nrho_N^AB(1N1N)\n=\nNPhi_A(1N)conjPhi_B(1N)\n=\nleftbeginmatrix\nrho_1(11)rho_1(1N)\nvdotsddotsvdots\nrho_1(N1)rho_1(NN)\nendmatrixright\nendaligned\nendequation"
},

{
    "location": "density_matrices/#References-1",
    "page": "Density matrices",
    "title": "References",
    "category": "section",
    "text": "Per-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. I. Physical Interpretations by Means of Density Matrices, Natural Spin-Orbitals, and Convergence Problems in the Method of Configurational Interaction. Physical Review, 97(6), 1474–1489. 10.1103/physrev.97.1474\nPer-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. II. Study of the Ordinary Hartree-Fock Approximation. Physical Review, 97(6), 1490–1508. 10.1103/physrev.97.1490\nPer-Olov Löwdin (1955). Quantum Theory of Many-Particle Systems. III. Extension of the Hartree-Fock Scheme to Include Degenerate Systems and Correlation Effects. Physical Review, 97(6), 1509–1520. 10.1103/physrev.97.1509"
},

{
    "location": "radial_orbitals/#",
    "page": "Radial orbitals",
    "title": "Radial orbitals",
    "category": "page",
    "text": ""
},

{
    "location": "radial_orbitals/#Atoms.RadialOrbital",
    "page": "Radial orbitals",
    "title": "Atoms.RadialOrbital",
    "category": "type",
    "text": "RadialOrbital\n\nA radial orbital is represented using basis coupled with a vector of expansion coefficients with respect to that basis; the basis implemented as an AbstractQuasimatrix.\n\n\n\n\n\n"
},

{
    "location": "radial_orbitals/#Atoms.RadialOrbitals",
    "page": "Radial orbitals",
    "title": "Atoms.RadialOrbitals",
    "category": "type",
    "text": "RadialOrbitals\n\nA collection of radial orbitals is instead represented using a matrix of such expansion coefficients, where each matrix column corresponds to a single radial orbital.\n\n\n\n\n\n"
},

{
    "location": "radial_orbitals/#Radial-orbitals-1",
    "page": "Radial orbitals",
    "title": "Radial orbitals",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendRadialOrbital\nRadialOrbitalsCurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "atom_types/#",
    "page": "Atom types",
    "title": "Atom types",
    "category": "page",
    "text": ""
},

{
    "location": "atom_types/#Atoms.Atom",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "type",
    "text": "Atom(radial_orbitals, orbitals, configurations, mix_coeffs, potential)\n\nAn atom constitutes a set of single-electron orbitals with associated radial_orbitals, configurations which are ManyElectronWavefunction:s, comprising of anti-symmetrized combinations of such orbitals. The expansion coefficients mix_coeffs determine the linear combination of the configurations for multi-configurational atoms.\n\nThe potential can be used to model either the nucleus by itself (a point charge or a nucleus of finite extent) or the core orbitals (i.e. a pseudo-potential).\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{O}, Tuple{B}, Tuple{T}, Tuple{MulQuasiArray{T,2,#s16} where #s16<:(Mul{#s15,#s14} where #s14<:(Tuple{#s13,#s12} where #s12<:(AbstractArray{T,2} where T) where #s13<:B) where #s15<:Tuple),Array{O,1},Array{#s1,1} where #s1<:TC,P,Type{C}}} where P where C where TC<:(Union{CSF{O,IT,T} where T<:Union{Term, HalfInteger} where IT<:Union{IntermediateTerm, HalfInteger}, Configuration{#s18} where #s18<:SpinOrbital} where O<:AtomicLevels.AbstractOrbital) where O where B where T<:Number",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(radial_orbitals, orbitals, configurations, potential, ::Type{C})\n\nCreate an Atom from the lists of radial_orbitals associated with orbitals and electronic configurations, with a nucleus modelled by potential. C determines the eltype of the mixing coefficients, which are initialized to [1,0,0,...].\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{B}, Tuple{T}, Tuple{UndefInitializer,Type{T},B,Array{TC,1},P,Type{C}}} where P where C where TC where B<:ContinuumArrays.QuasiArrays.AbstractQuasiArray{T,2} where T<:Number",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(undef, ::Type{T}, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})\n\nCreate an Atom on the space spanned by R, from the list of electronic configurations, with a nucleus modelled by potential, and leave the orbitals uninitialized. T determines the eltype of the radial orbitals and C the mixing coefficients.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{B}, Tuple{T}, Tuple{Symbol,Type{T},B,Array{TC,1},P,Type{C}}} where P where C where TC where B<:ContinuumArrays.QuasiArrays.AbstractQuasiArray{T,2} where T<:Number",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(init, ::Type{T}, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})\n\nCreate an Atom on the space spanned by R, from the list of electronic configurations, with a nucleus modelled by potential, and initialize the orbitals according to init. T determines the eltype of the radial orbitals and C the mixing coefficients.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{B}, Tuple{T}, Tuple{Init}, Tuple{Init,B,Array{TC,1},P,Type{C}}} where P<:AtomicPotentials.AbstractPotential where C where TC where B<:ContinuumArrays.QuasiArrays.AbstractQuasiArray{T,2} where T where Init",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(init, R::AbstractQuasiMatrix, configurations, potential, ::Type{C})\n\nCreate an Atom on the space spanned by R, from the list of electronic configurations, with a nucleus modelled by potential, and initialize the orbitals according to init. C determines the eltype of the mixing coefficients.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{B}, Tuple{B,Array{#s1,1} where #s1<:TC,P}, Tuple{B,Array{#s2,1} where #s2<:TC,P,Type{C}}} where P<:AtomicPotentials.AbstractPotential where C where TC<:(Union{CSF{O,IT,T} where T<:Union{Term, HalfInteger} where IT<:Union{IntermediateTerm, HalfInteger}, Configuration{#s18} where #s18<:SpinOrbital} where O<:AtomicLevels.AbstractOrbital) where B<:(ContinuumArrays.QuasiArrays.AbstractQuasiArray{T,2} where T)",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(R::AbstractQuasiMatrix, configurations, potential[, ::Type{C}=eltype(R)])\n\nCreate an Atom on the space spanned by R, from the list of electronic configurations, with a nucleus modelled by potential, and initialize the orbitals to their hydrogenic values.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.Atom-Union{Tuple{P}, Tuple{C}, Tuple{TC}, Tuple{O}, Tuple{B}, Tuple{T}, Tuple{Atom{T,B,O,TC,C,P},Array{#s2,1} where #s2<:TC}} where P where C where TC where O where B where T",
    "page": "Atom types",
    "title": "Atoms.Atom",
    "category": "method",
    "text": "Atom(other_atom::Atom, configurations)\n\nCreate a new atom using the same basis and nuclear potential as other_atom, but with a different set of configurations. The orbitals of other_atom are copied over as starting guess.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.DiracAtom",
    "page": "Atom types",
    "title": "Atoms.DiracAtom",
    "category": "type",
    "text": "DiracAtom\n\nA DiracAtom is a specialization of Atom for the relativistic case.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Base.getindex",
    "page": "Atom types",
    "title": "Base.getindex",
    "category": "function",
    "text": "getindex(atom, j)\n\nReturns a copy of the j:th radial orbital.\n\n\n\n\n\ngetindex(atom, orb)\n\nReturns a copy of the radial orbital corresponding to orb.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Base.view",
    "page": "Atom types",
    "title": "Base.view",
    "category": "function",
    "text": "view(atom, j)\n\nReturns a view of the j:th radial orbital.\n\n\n\n\n\nview(atom, orb)\n\nReturns a view of the radial orbital corresponding to orb.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#AtomicLevels.num_electrons",
    "page": "Atom types",
    "title": "AtomicLevels.num_electrons",
    "category": "function",
    "text": "num_electrons(atom)\n\nReturn number of electrons in atom.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atom-types-1",
    "page": "Atom types",
    "title": "Atom types",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendAtom\nAtom(radial_orbitals::RadialOrbitals{T,B}, orbitals::Vector{O}, configurations::Vector{<:TC}, potential::P, ::Type{C}) where {T<:Number,B,O,TC<:ManyElectronWavefunction,C,P}\nAtom(::UndefInitializer, ::Type{T}, R::B, configurations::Vector{TC}, potential::P, ::Type{C}) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}\nAtom(init::Symbol, ::Type{T}, R::B, configurations::Vector{TC}, potential::P, ::Type{C}; kwargs...) where {T<:Number,B<:AbstractQuasiMatrix{T},TC,C,P}\nAtom(init::Init, R::B, configurations::Vector{TC}, potential::P, ::Type{C}; kwargs...) where {Init,T,B<:AbstractQuasiMatrix{T},TC,C,P<:AbstractPotential}\nAtom(R::B, configurations::Vector{<:TC}, potential::P, ::Type{C}=eltype(R); kwargs...) where {B<:AbstractQuasiMatrix,TC<:ManyElectronWavefunction,C,P<:AbstractPotential}\nAtom(other_atom::Atom{T,B,O,TC,C,P}, configurations::Vector{<:TC}; kwargs...) where {T,B,O,TC,C,P}\nDiracAtomgetindex\nview\nnum_electrons"
},

{
    "location": "atom_types/#Atoms.ManyElectronWavefunction",
    "page": "Atom types",
    "title": "Atoms.ManyElectronWavefunction",
    "category": "constant",
    "text": "ManyElectronWavefunction\n\nA many-electron wave function configuration can either be given as a CSF (summed over spins) or a configuration of spin-orbitals (where all quantum numbers are specified).\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.outsidecoremodel",
    "page": "Atom types",
    "title": "Atoms.outsidecoremodel",
    "category": "function",
    "text": "outsidecoremodel(configuration::Configuration, potential::P)\n\nReturn the part of the electronic configuration that is not part of the the configuration modelled by the potential. For a point charge, this is the same as the configuration itself, but for pseudopotential, typically only the outer shells remain.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Atoms.all_bound",
    "page": "Atom types",
    "title": "Atoms.all_bound",
    "category": "function",
    "text": "all_bound(atom)\n\nReturns true if all orbitals in atom are bound orbitals.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#SCF.coefficients",
    "page": "Atom types",
    "title": "SCF.coefficients",
    "category": "function",
    "text": "SCF.coefficients(atom)\n\nReturns a view of the mixing coefficients.\n\n\n\n\n\n"
},

{
    "location": "atom_types/#SCF.orbitals",
    "page": "Atom types",
    "title": "SCF.orbitals",
    "category": "function",
    "text": "SCF.orbitals(atom)\n\nReturns a view of the radial orbital coefficients (NB, it does not return the MulQuasiMatrix, but the actual underlying expansion coefficients, since SCF operates on them in the self-consistent iteration).\n\n\n\n\n\n"
},

{
    "location": "atom_types/#Internals-1",
    "page": "Atom types",
    "title": "Internals",
    "category": "section",
    "text": "CurrentModule = AtomsManyElectronWavefunction\noutsidecoremodel\nall_bound\nSCF.coefficients\nSCF.orbitalsCurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "one_body/#",
    "page": "One-body Hamiltonians",
    "title": "One-body Hamiltonians",
    "category": "page",
    "text": ""
},

{
    "location": "one_body/#Atoms.one_body_hamiltonian",
    "page": "One-body Hamiltonians",
    "title": "Atoms.one_body_hamiltonian",
    "category": "function",
    "text": "one_body_hamiltonian(::Type{Tuple}, atom, orb)\n\nReturn the kinetic and one-body potential energy operators (as a tuple) for the orbital orb of atom.\n\n\n\n\n\none_body_hamiltonian(::Type{Tuple}, atom, orb)\n\nReturn the one-body energy operator for the orbital orb of atom.\n\n\n\n\n\n"
},

{
    "location": "one_body/#Atoms.KineticEnergyHamiltonian",
    "page": "One-body Hamiltonians",
    "title": "Atoms.KineticEnergyHamiltonian",
    "category": "type",
    "text": "KineticEnergyHamiltonian\n\nThe kinetic energy part of the one-body Hamiltonian, excluding the centrifugal potential. It is diagonal in spin, i.e. it does not couple orbitals of opposite spin.\n\n\n\n\n\n"
},

{
    "location": "one_body/#Atoms.PotentialEnergyHamiltonian",
    "page": "One-body Hamiltonians",
    "title": "Atoms.PotentialEnergyHamiltonian",
    "category": "type",
    "text": "PotentialEnergyHamiltonian\n\nThe potential energy part of the one-body Hamiltonian, including the centrifugal potential. It is diagonal in spin, i.e. it does not couple orbitals of opposite spin.\n\n\n\n\n\n"
},

{
    "location": "one_body/#Atoms.AtomicOneBodyHamiltonian",
    "page": "One-body Hamiltonians",
    "title": "Atoms.AtomicOneBodyHamiltonian",
    "category": "type",
    "text": "AtomicOneBodyHamiltonian(op, orbital)\n\nStructure holding a one-body energy operator op acting on its associated orbital.\n\n\n\n\n\n"
},

{
    "location": "one_body/#LazyArrays.:⋆",
    "page": "One-body Hamiltonians",
    "title": "LazyArrays.:⋆",
    "category": "function",
    "text": "ĥ ⋆ ϕ\n\nReturn the lazy multiplicative action of the AtomicOneBodyHamiltonian ĥ on the radial orbital coefficient vector ϕ.\n\n\n\n\n\n"
},

{
    "location": "one_body/#Base.Broadcast.materialize!-Union{Tuple{MulAdd{#s4,#s3,#s2,T,#s1,Source,Dest} where #s1<:AtomicOneBodyHamiltonian where #s2 where #s3 where #s4}, Tuple{Dest}, Tuple{Source}, Tuple{T}} where Dest where Source where T",
    "page": "One-body Hamiltonians",
    "title": "Base.Broadcast.materialize!",
    "category": "method",
    "text": "materialize!(::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*H*x + β*y where H is a AtomicOneBodyHamiltonian and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a DirectPotential (with a precomputed direct potential computed via SCF.update!) and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a ExchangePotential (by solving the Poisson problem with x as one of the constituent source orbitals in the mutual density) and x and y are RadialOrbitals.\n\n\n\n\n\n"
},

{
    "location": "one_body/#One-body-Hamiltonians-1",
    "page": "One-body Hamiltonians",
    "title": "One-body Hamiltonians",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendThe one-body Hamiltonian for electron i in an atom is given bybeginequation\nhamiltonian_i defd\n-fracnabla_i^22 +\nV(vecr_i)\nendequationIn spherical coordinates (and using reduced wavefunctions), the Laplacian transforms tobeginequation\noperatorT_r_i = -fracpartial_r_i^22 + fracell(ell+1)2r_i^2\nendequationwhere the second term, called the centrifugal potential, although originating from the Laplacian, is usually treated together with the nuclear potential V(r_i).one_body_hamiltonian\nKineticEnergyHamiltonian\nPotentialEnergyHamiltonian\nAtomicOneBodyHamiltonian\nLazyArrays.:(⋆)\nLazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest}) where {T,Source,Dest}"
},

{
    "location": "one_body/#Atoms.diagonalize_one_body",
    "page": "One-body Hamiltonians",
    "title": "Atoms.diagonalize_one_body",
    "category": "function",
    "text": "diagonalize_one_body(H, nev; method=:arnoldi_shift_invert, tol=1e-10, σ=-1)\n\nDiagonalize the one-body Hamiltonian H and find the nev lowest eigenpairs, using the specified diagonalization method; valid choices are\n\n:arnoldi which performs the standard Krylov iteration looking for the eigenvalues with smallest real values,\n:arnoldi_shift_invert which performs the Krylov iteration but with the shifted and inverted matrix (H - I*σ)⁻¹ looking for the eigenvalues with largest real values,\n:eigen which uses Julia\'s built-in eigensolver.\n\ntol sets the Krylov tolerance.\n\n\n\n\n\n"
},

{
    "location": "one_body/#Diagonalization-of-one-body-Hamiltonians-1",
    "page": "One-body Hamiltonians",
    "title": "Diagonalization of one-body Hamiltonians",
    "category": "section",
    "text": "diagonalize_one_bodyCurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "hydrogenic/#",
    "page": "Hydrogenic initialization",
    "title": "Hydrogenic initialization",
    "category": "page",
    "text": ""
},

{
    "location": "hydrogenic/#Atoms.hydrogenic!",
    "page": "Hydrogenic initialization",
    "title": "Atoms.hydrogenic!",
    "category": "function",
    "text": "hydrogenic!(atom[; kwargs...])\n\nInitialize the radial orbitals of atom to their unscreened hydrogenic values. This is done via simple diagonalization of the one-body Hamiltonian for each angular symmetry. The kwargs are passed on to diagonalize_one_body and can be used to influence how the diagonalization is performed.\n\n\n\n\n\n"
},

{
    "location": "hydrogenic/#Hydrogenic-initialization-1",
    "page": "Hydrogenic initialization",
    "title": "Hydrogenic initialization",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendhydrogenic!CurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "orbital_equations/#",
    "page": "Orbital equations",
    "title": "Orbital equations",
    "category": "page",
    "text": ""
},

{
    "location": "orbital_equations/#Orbital-equations-1",
    "page": "Orbital equations",
    "title": "Orbital equations",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nend"
},

{
    "location": "orbital_equations/#Atoms.OrbitalHamiltonianTerm",
    "page": "Orbital equations",
    "title": "Atoms.OrbitalHamiltonianTerm",
    "category": "type",
    "text": "OrbitalHamiltonianTerm(i, j, coeff, A, integrals)\n\nRepresents a term in the orbital Hamiltonian arising from a variation of the energy expressions between configurations i and j in the multi-configurational expansion. coeff is the numeric coefficient, A is the operator acting on the orbital, and integrals is a vector of OrbitalIntegrals arising from the presence of non-orthogonal orbitals and whose values should be multiplied to form the overall coefficient.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.coefficient",
    "page": "Orbital equations",
    "title": "Atoms.coefficient",
    "category": "function",
    "text": "coefficient(term::OrbitalHamiltonianTerm)\n\nReturn the multiplicative coefficient pertaining to term, excluding the conj(c_i)*c_j mixing coefficients, due to the configuration-interaction.\n\n\n\n\n\ncoefficient(term::OrbitalHamiltonianTerm, c::Vector)\n\nReturn the multiplicative coefficient pertaining to term, including the conj(c_i)*c_j mixing coefficients, due to the configuration-interaction.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.OrbitalHamiltonian",
    "page": "Orbital equations",
    "title": "Atoms.OrbitalHamiltonian",
    "category": "type",
    "text": "OrbitalHamiltonian(R, terms, mix_coeffs, projector, orbital)\n\nThe Hamiltonian for orbital is constructed from a radial basis R, a set of OrbitalHamiltonianTerm terms that describe the various interactions between orbitals, mix_coeffs which are the mixing coefficents for the multi-configurational expansion. The projector ensures orthogonality between orbital pairs which have Lagrange multipliers associated with them, by projecting out components of other orbitals every time the OrbitalHamiltonian action on orbital is computed.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.Projector",
    "page": "Orbital equations",
    "title": "Atoms.Projector",
    "category": "type",
    "text": "Projector(ϕs)\n\nRepresents the projector out of the subspace spanned by the radial orbitals ϕs\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.projectout!",
    "page": "Orbital equations",
    "title": "Atoms.projectout!",
    "category": "function",
    "text": "projectout!(y, projector)\n\nProject out all components of y parallel to the radial orbitals projector.ϕs.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#SCF.energy_matrix!-Union{Tuple{B}, Tuple{T}, Tuple{O}, Tuple{HM}, Tuple{HM,OrbitalHamiltonian{O,T,B,OV,Proj} where Proj where OV,MulQuasiArray{T,1,#s12} where #s12<:(Mul{#s13,#s14} where #s14<:(Tuple{#s15,#s16} where #s16<:(AbstractArray{T,1} where T) where #s15<:B) where #s13<:Tuple)}} where B where T where O where HM<:(AbstractArray{T,2} where T)",
    "page": "Orbital equations",
    "title": "SCF.energy_matrix!",
    "category": "method",
    "text": "energy_matrix!(H, hamiltonian, ϕ)\n\nCompute the contribution of hamiltonian to the Hamiltonian matrix H by repeatedly acting on the associated radial orbital ϕ with the different multi-configurational OrbitalHamiltonianTerms of hamiltonian.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.filter-Tuple{Function,Atoms.OrbitalHamiltonian}",
    "page": "Orbital equations",
    "title": "Base.filter",
    "category": "method",
    "text": "filter(fun::Function, H::OrbitalHamiltonian)\n\nFilter the OrbitalHamiltonianTerms of H according to the predicate fun.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.copyto!-Union{Tuple{M}, Tuple{T}, Tuple{M,OrbitalHamiltonian}} where M<:AbstractArray{T,2} where T",
    "page": "Orbital equations",
    "title": "Base.copyto!",
    "category": "method",
    "text": "copyto!(dest::AbstractMatix, hamiltonian::OrbitalHamiltonian)\n\nMaterialize the orbital hamiltonian into matrix form and store it in dest, using the current values of all other orbitals. This is only possible if the orbital hamiltonian does not contain any ExchangePotentials or SourceTerms, since the former is non-local (and thus not representable as a matrix) and the latter is not a linear operator (but an affine one).\n\nTypical usage is to compute an easily factorizable matrix that can be used for preconditioning the solution of the full equation.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.:+-Union{Tuple{Proj}, Tuple{OV}, Tuple{B}, Tuple{T}, Tuple{O}, Tuple{OrbitalHamiltonian{O,T,B,OV,Proj},UniformScaling}} where Proj where OV where B where T where O",
    "page": "Orbital equations",
    "title": "Base.:+",
    "category": "method",
    "text": "h::OrbitalHamiltonian + λ::UniformScaling\n\nShift the OrbitalHamiltonian h by λ.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.:--Tuple{Atoms.OrbitalHamiltonian,LinearAlgebra.UniformScaling}",
    "page": "Orbital equations",
    "title": "Base.:-",
    "category": "method",
    "text": "h::OrbitalHamiltonian - λ::UniformScaling\n\nShift the OrbitalHamiltonian h by -λ.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#SCF.KrylovWrapper",
    "page": "Orbital equations",
    "title": "SCF.KrylovWrapper",
    "category": "type",
    "text": "SCF.KrylovWrapper(hamiltonian::OrbitalHamiltonian)\n\nConstruct a KrylovWrapper such that hamiltonian, that acts on function spaces, can be used in a Krylov solver, which works with linear algebra vector spaces.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#LinearAlgebra.mul!-Union{Tuple{Hamiltonian}, Tuple{B}, Tuple{T}, Tuple{V₂}, Tuple{V₁}, Tuple{V₁,KrylovWrapper{T,Hamiltonian},V₂}} where Hamiltonian<:Atoms.OrbitalHamiltonian where B where T where V₂ where V₁",
    "page": "Orbital equations",
    "title": "LinearAlgebra.mul!",
    "category": "method",
    "text": "mul!(y, A::KrylovWrapper{T,<:OrbitalHamiltonian}, x)\n\nMaterialize the action of the OrbitalHamiltonian on the linear algebra vector x and store the result in y, by wrapping them both with the QuasiMatrix necessary to transform x and y to the function space of the Hamiltonian.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#MatrixFactorizations.preconditioner",
    "page": "Orbital equations",
    "title": "MatrixFactorizations.preconditioner",
    "category": "function",
    "text": "MatrixFactorizations.preconditioner(hamiltonian::OrbitalHamiltonian)\n\nReturn a factorization of the matrix corresponding to hamiltonian, where all terms arising from exchange and configuration interaction have been removes, since they cannot be represented by a matrix.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Hamiltonian-1",
    "page": "Orbital equations",
    "title": "Hamiltonian",
    "category": "section",
    "text": "OrbitalHamiltonianTerm\ncoefficient\nOrbitalHamiltonian\nProjector\nprojectout!\nSCF.energy_matrix!(H::HM, hamiltonian::OrbitalHamiltonian{O,T,B}, ϕ::RadialOrbital{T,B}) where {HM<:AbstractMatrix,O,T,B}\nBase.filter(fun::Function, H::OrbitalHamiltonian)\nBase.copyto!(dest::M, hamiltonian::OrbitalHamiltonian) where {T,M<:AbstractMatrix{T}}\nBase.:(+)(h::OrbitalHamiltonian{O,T,B,OV,Proj}, λ::UniformScaling) where {O,T,B,OV,Proj}\nBase.:(-)(h::OrbitalHamiltonian, λ::UniformScaling)\nSCF.KrylovWrapper\nLinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,Hamiltonian}, x::V₂) where {V₁,V₂,T,B,Hamiltonian<:OrbitalHamiltonian}\nMatrixFactorizations.preconditioner"
},

{
    "location": "orbital_equations/#Atoms.OrbitalIntegral",
    "page": "Orbital equations",
    "title": "Atoms.OrbitalIntegral",
    "category": "type",
    "text": "OrbitalIntegral{N}\n\nAbstract type for integrals of rank N of orbitals, whose values need to be recomputed every time the orbitals are updated. Rank 0 corresponds to a scalar value, rank 1 to a diagonal matrix, etc.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.OrbitalOverlapIntegral",
    "page": "Orbital equations",
    "title": "Atoms.OrbitalOverlapIntegral",
    "category": "type",
    "text": "OrbitalOverlapIntegral(a, b, av, bv, value)\n\nRepresents the orbital overlap integral ⟨a|b⟩, for orbitals a and b, along with views of their radial orbitals av and bv and the current value of the integral.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#SCF.update!-Tuple{Atoms.OrbitalOverlapIntegral}",
    "page": "Orbital equations",
    "title": "SCF.update!",
    "category": "method",
    "text": "SCF.update!(oo::OrbitalOverlapIntegral)\n\nUpdate the value of the integral oo.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.HFPotential",
    "page": "Orbital equations",
    "title": "Atoms.HFPotential",
    "category": "type",
    "text": "HFPotential(k, a, b, av, bv, V̂, poisson)\n\nRepresents the k:th multipole exansion of the Hartree–Fock potential formed by orbitals a and b (av and bv being views of their corresponding radial orbitals). V̂ is the resultant one-body potential formed, which can act on a third orbital and poisson computes the potential by solving Poisson\'s problem.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.DirectPotential",
    "page": "Orbital equations",
    "title": "Atoms.DirectPotential",
    "category": "type",
    "text": "DirectPotential\n\nSpecial case of HFPotential for the direct interaction, in which case the potential formed from two orbitals can be precomputed before acting on a third orbital.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#SCF.update!-Union{Tuple{HFPotential{:direct,O,T,B,OV,RO,P}}, Tuple{P}, Tuple{RO}, Tuple{OV}, Tuple{B}, Tuple{T}, Tuple{O}} where P where RO where OV where B where T where O",
    "page": "Orbital equations",
    "title": "SCF.update!",
    "category": "method",
    "text": "SCF.update!(p::DirectPotential)\n\nUpdate the direct potential p by solving the Poisson problem with the current values of the orbitals forming the mutual density.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.Broadcast.materialize!-Union{Tuple{MulAdd{#s2,#s1,#s3,T,#s4,Source,Dest} where #s4<:(HFPotential{:direct,O,T,B,OV,RO,P} where P where RO where OV where B where T where O) where #s3 where #s1 where #s2}, Tuple{Dest}, Tuple{Source}, Tuple{T}} where Dest where Source where T",
    "page": "Orbital equations",
    "title": "Base.Broadcast.materialize!",
    "category": "method",
    "text": "materialize!(::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*H*x + β*y where H is a AtomicOneBodyHamiltonian and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a DirectPotential (with a precomputed direct potential computed via SCF.update!) and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a ExchangePotential (by solving the Poisson problem with x as one of the constituent source orbitals in the mutual density) and x and y are RadialOrbitals.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.ExchangePotential",
    "page": "Orbital equations",
    "title": "Atoms.ExchangePotential",
    "category": "type",
    "text": "ExchangePotential\n\nSpecial case of HFPotential for the exchange interaction, in which case the potential is formed from the orbital acted upon, along with another orbital, and then applied to a third orbital. Thus this potential cannot be precomputed, but must be recomputed every time the operator is applied. This makes this potential expensive to handle and the number of times it is applied should be minimized, if possible.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Base.Broadcast.materialize!-Union{Tuple{MulAdd{#s4,#s3,#s2,T,#s1,Source,Dest} where #s1<:(HFPotential{:exchange,O,T,B,OV,RO,P} where P where RO where OV where B where T where O) where #s2 where #s3 where #s4}, Tuple{Dest}, Tuple{Source}, Tuple{T}} where Dest where Source where T",
    "page": "Orbital equations",
    "title": "Base.Broadcast.materialize!",
    "category": "method",
    "text": "materialize!(::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*H*x + β*y where H is a AtomicOneBodyHamiltonian and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a DirectPotential (with a precomputed direct potential computed via SCF.update!) and x and y are RadialOrbitals.\n\n\n\n\n\nmaterialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*V̂*x + β*y where V̂ is a ExchangePotential (by solving the Poisson problem with x as one of the constituent source orbitals in the mutual density) and x and y are RadialOrbitals.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.SourceTerm",
    "page": "Orbital equations",
    "title": "Atoms.SourceTerm",
    "category": "type",
    "text": "SourceTerm(operator, source_orbital, ov)\n\nThe point of SourceTerm is to implement inhomogeneous terms that contribute to the equation for an orbital, and whose input is some other source_orbital. This kind of term appears in multi-configurational problems.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Atoms.ShiftTerm",
    "page": "Orbital equations",
    "title": "Atoms.ShiftTerm",
    "category": "type",
    "text": "ShiftTerm(λ)\n\nThe point of ShiftTerm is to implement an overall energy shift of the Hamiltonian.\n\n\n\n\n\n"
},

{
    "location": "orbital_equations/#Orbital-integrals-and-terms-1",
    "page": "Orbital equations",
    "title": "Orbital integrals and terms",
    "category": "section",
    "text": "OrbitalIntegral\nOrbitalOverlapIntegral\nSCF.update!(oo::OrbitalOverlapIntegral; kwargs...)\nHFPotential\nDirectPotential\nSCF.update!(p::DirectPotential{O,T,B,OV,RO,P}; kwargs...) where {O,T,B,OV,RO,P}\nLazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:DirectPotential, Source, Dest}) where {T,Source,Dest}\nExchangePotential\nLazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:ExchangePotential, Source, Dest}) where {T,Source,Dest}\nSourceTerm\nShiftTerm"
},

{
    "location": "orbital_equations/#Orbital-equations-2",
    "page": "Orbital equations",
    "title": "Orbital equations",
    "category": "section",
    "text": "AtomicOrbitalEquation\nenergy(hfeq::AtomicOrbitalEquation, term=:all)CurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "equation_systems/#",
    "page": "Equation systems",
    "title": "Equation systems",
    "category": "page",
    "text": ""
},

{
    "location": "equation_systems/#Atoms.AtomicEquations",
    "page": "Equation systems",
    "title": "Atoms.AtomicEquations",
    "category": "type",
    "text": "AtomicEquations(atom, equations, integrals)\n\nStructure representing the (e.g. Hartree–Fock) equations for atom, along with all integrals that are shared between the equations.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#SCF.update!-Tuple{Atoms.AtomicEquations}",
    "page": "Equation systems",
    "title": "SCF.update!",
    "category": "method",
    "text": "update!(equations::AtomicEquations)\n\nRecompute all integrals using the current values for the radial orbitals.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#SCF.energy_matrix!-Union{Tuple{HM}, Tuple{HM,AtomicEquations}, Tuple{HM,AtomicEquations,Symbol}} where HM<:(AbstractArray{T,2} where T)",
    "page": "Equation systems",
    "title": "SCF.energy_matrix!",
    "category": "method",
    "text": "energy_matrix!(H, hfeqs::AtomicEquations[, which=:energy])\n\nCompute the energy matrix by computing the energy observable and storing it in H. Requires that hfeqs has the :energy and :kinetic_energy Observables registered (this is the default).\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Atoms.find_symmetries",
    "page": "Equation systems",
    "title": "Atoms.find_symmetries",
    "category": "function",
    "text": "find_symmetries(orbitals)\n\nGroup all orbitals according to their symmetries, e.g. ℓ for Orbitals. This is used to determine which off-diagonal Lagrange multipliers are necessary to maintain orthogonality.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Atoms.generate_atomic_orbital_equations",
    "page": "Equation systems",
    "title": "Atoms.generate_atomic_orbital_equations",
    "category": "function",
    "text": "generate_atomic_orbital_equations(atom::Atom, eqs::MCEquationSystem,\n                                  integrals, integral_map)\n\nFor each variationally derived orbital equation in eqs, generate the corresponding AtomicOrbitalEquation.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Base.diff",
    "page": "Equation systems",
    "title": "Base.diff",
    "category": "function",
    "text": "diff(atom[, H]; overlaps=[], selector=outsidecoremodel, verbosity=0)\n\nDifferentiate the energy expression of the Hamiltonian H associated with the atom\'s configurations(s) with respect to the atomic orbitals to derive the Hartree–Fock equations for the orbitals.\n\nBy default, the Hamiltonian H=FieldFreeOneBodyHamiltonian()+CoulombInteraction().\n\nNon-orthogonality between orbitals can be specified by providing OrbitalOverlaps between these pairs. Only those electrons not modelled by atom.potential of each configuration are considered for generating the energy expression, this can be changed by choosing another value for selector.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Equation-systems-1",
    "page": "Equation systems",
    "title": "Equation systems",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendAtomicEquations\nSCF.update!(equations::AtomicEquations; kwargs...)\nSCF.energy_matrix!(H::HM, hfeqs::AtomicEquations, which::Symbol=:total) where {HM<:AbstractMatrix}\nfind_symmetries\ngenerate_atomic_orbital_equations\nBase.diff"
},

{
    "location": "equation_systems/#Atoms.push_common_integral!",
    "page": "Equation systems",
    "title": "Atoms.push_common_integral!",
    "category": "function",
    "text": "push_common_integral!(integrals, integral_map,\n                      integral, atom)\n\nPush the integral to integrals, constructing the correct OrbitalIntegral pertaining to the atom. Record the index of integral in integrals in the integral_map.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Atoms.pushterms!",
    "page": "Equation systems",
    "title": "Atoms.pushterms!",
    "category": "function",
    "text": "pushterms!(terms, operator, equation_terms,\n           integrals, integral_map, symbolic_integrals)\n\nFor each term in equation_terms, push a term, located at CI coordinates i,j, of the overall orbital Hamiltonian to terms, constructed from operator and a product of orbital integrals, multiplied by an overall factor given by expression and multipole expansions. integrals contain common OrbitalIntegrals and integral_map maps from symbolic_integrals to integrals.\n\n\n\n\n\n"
},

{
    "location": "equation_systems/#Common-integrals-1",
    "page": "Equation systems",
    "title": "Common integrals",
    "category": "section",
    "text": "When deriving the equations of motion from an energy expression, the same integral may appear many times in the equations for different orbitals, multiplied by different factors, etc. To minimize the reevaluation of integrals, AtomicEquations keeps track of all the common integrals, and they are recomputed exactly once, when SCF.update! is called. The routines below are used when setting up the equation system.push_common_integral!\npushterms!CurrentModule = nothing\nDocTestSetup = nothing"
},

{
    "location": "observables/#",
    "page": "Observables",
    "title": "Observables",
    "category": "page",
    "text": ""
},

{
    "location": "observables/#Atoms.Observable",
    "page": "Observables",
    "title": "Atoms.Observable",
    "category": "type",
    "text": "Observable\n\nRepresents a physical quantity that can be observed, which is calculated as the matrix element of an operator between two configurations. All physical observables are real.\n\n\n\n\n\n"
},

{
    "location": "observables/#Atoms.Observable-Union{Tuple{A}, Tuple{QuantumOperator,A,Array{#s4,1} where #s4<:OrbitalOverlap,Array{OrbitalIntegral,1},Dict{Any,Int64},Dict}} where A<:Atom",
    "page": "Observables",
    "title": "Atoms.Observable",
    "category": "method",
    "text": "Observable(operator, atom, overlaps, integrals)\n\nConstruct an observable corresponding the operator acting on atom. overlaps is a list of non-orthogonal, integrals a list of common integrals, and integral_map is a mapping from symbolic integrals to OrbitalIntegrals.\n\n\n\n\n\n"
},

{
    "location": "observables/#Atoms.observe!",
    "page": "Observables",
    "title": "Atoms.observe!",
    "category": "function",
    "text": "observe!(A::M, o::Observable)\n\nCompute the observable o between all configurations and store the results as matrix elements of A.\n\n\n\n\n\n"
},

{
    "location": "observables/#Observables-1",
    "page": "Observables",
    "title": "Observables",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendObservable\nObservable(operator::QuantumOperator, atom::A, overlaps::Vector{<:OrbitalOverlap}, integrals::Vector{OrbitalIntegral}, integral_map::Dict{Any,Int}, symmetries::Dict) where {A<:Atom}\nobserve!CurrentModule = nothing\nDocTestSetup = nothing"
},

]}
