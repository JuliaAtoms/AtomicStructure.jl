module SCF

using LinearAlgebra
using BandedMatrices
using SparseArrays
using IntervalSets

using ArnoldiMethod
using IterativeSolvers
using Optim
using LineSearches

using SolverTraces
import SolverTraces: base_exp

using UnicodeFun
using Formatting
using Crayons

using Compat

include("quantum_systems.jl")
include("krylov_wrapper.jl")
include("fock.jl")
include("solver_trace.jl")
include("utils.jl")
include("secular_problem.jl")
include("gershgorin_discs.jl")
include("rotate.jl")
include("solve_orbital_equation.jl")
include("self_consistent_iteration.jl")
include("fock_problem.jl")
include("manifolds.jl")
include("optim.jl")

export Fock, scf!

end # module
