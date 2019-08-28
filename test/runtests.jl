using Atoms
using CoulombIntegrals
import CoulombIntegrals: locs
using AtomicLevels
using HalfIntegers
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

function get_atom_grid(grid_type, rₘₐₓ, ρ, nucleus; fedvr_order=10, amend_order=true)
    Z = charge(nucleus)

    R,r = if grid_type == :fedvr
        N = max(ceil(Int, rₘₐₓ/(ρ*fedvr_order)),2)
        t = range(0.0, stop=rₘₐₓ, length=N)
        amended_order = vcat(fedvr_order+5, Fill(fedvr_order,length(t)-2))
        FEDVR(t, amend_order ? amended_order : fedvr_order)[:,2:end-1], range(t[1],stop=t[end],length=1001)
    else
        N = ceil(Int, rₘₐₓ/ρ + 1/2)
        R=RadialDifferences(N, ρ, Z)
        R,FiniteDifferencesQuasi.locs(R)
    end

    R,r
end

include("test_slater_integrals.jl")
include("structure_setup.jl")
include("calculation_accuracy.jl")
include("spin_orbit.jl")
