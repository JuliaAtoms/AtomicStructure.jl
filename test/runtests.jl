using Atoms
using AtomicLevels
using AtomicPotentials
using PseudoPotentials
using AngularMomentumAlgebra
using SCF
using PrettyTables

# Needed as /an/ implementation of an AbstractQuasiMatrix, not the
# only possibility:
using FiniteDifferencesQuasi
using FEDVRQuasi

using FillArrays

using Test

include("test_slater_integrals.jl")
include("structure_setup.jl")
include("calculation_accuracy.jl")
