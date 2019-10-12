using Documenter
using EnergyExpressions
using AtomicLevels
using Atoms

DocMeta.setdocmeta!(Atoms, :DocTestSetup, :(using Atoms, AtomicLevels); recursive=true)
makedocs(
    modules = [Atoms],
    sitename = "Atoms",
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Density matrices" => "density_matrices.md"
        ],
        "Implementation" => [
            "Radial orbitals" => "radial_orbitals.md",
            "Atom types" => "atom_types.md",
            "One-body Hamiltonians" => "one_body.md",
            "Hydrogenic initialization" => "hydrogenic.md",
            "Orbital equations" => "orbital_equations.md",
            "Equation systems" => "equation_systems.md",
            "Observables" => "observables.md"
        ]
    ],
    format = Documenter.HTML(assets = ["assets/latex.js"])
)

deploydocs(repo = "github.com/JuliaAtoms/Atoms.jl.git")
