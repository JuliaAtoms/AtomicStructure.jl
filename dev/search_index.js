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
    "location": "one_body/#Base.Broadcast.materialize!-Union{Tuple{MulAdd{#s2,#s1,#s3,T,#s4,Source,Dest} where #s4<:AtomicOneBodyHamiltonian where #s3 where #s1 where #s2}, Tuple{Dest}, Tuple{Source}, Tuple{T}} where Dest where Source where T",
    "page": "One-body Hamiltonians",
    "title": "Base.Broadcast.materialize!",
    "category": "method",
    "text": "materialize!(::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest})\n\nMaterialize the lazy multiplication–addition of the type y ← α*H*x + β*y where H is a AtomicOneBodyHamiltonian and x and y are RadialOrbitals.\n\n\n\n\n\n"
},

{
    "location": "one_body/#One-body-Hamiltonians-1",
    "page": "One-body Hamiltonians",
    "title": "One-body Hamiltonians",
    "category": "section",
    "text": "CurrentModule = Atoms\nDocTestSetup = quote\n    using Atoms\nendThe one-body Hamiltonian for electron i in an atom is given bybeginequation\nhamiltonian_i defd\n-fracnabla_i^22 +\nV(vecr_i)\nendequationIn spherical coordinates (and using reduced wavefunctions), the Laplacian transforms tobeginequation\noperatorT_r_i = -fracpartial_r_i^22 + fracell(ell+1)2r_i^2\nendequationwhere the second term, called the centrifugal potential, although originating from the Laplacian, is usually treated together with the nuclear potential V(r_i).one_body_hamiltonian\nKineticEnergyHamiltonian\nPotentialEnergyHamiltonian\nAtomicOneBodyHamiltonian\nLazyArrays.:(⋆)\nLazyArrays.materialize!(ma::MulAdd{<:Any, <:Any, <:Any, T, <:AtomicOneBodyHamiltonian, Source, Dest}) where {T,Source,Dest}"
},

{
    "location": "one_body/#Atoms.ShiftInvert",
    "page": "One-body Hamiltonians",
    "title": "Atoms.ShiftInvert",
    "category": "type",
    "text": "ShiftInvert(A⁻¹)\n\nHelp structure used in diagonalization of A via the shift-and-invert Krylov technique, where the action of A⁻¹ instead of A is computed in the Krylov iterations. This is useful for converging interior eigenvalues.\n\n\n\n\n\n"
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
    "text": "ShiftInvert\ndiagonalize_one_bodyCurrentModule = nothing\nDocTestSetup = nothing"
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

]}
