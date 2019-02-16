using Documenter
using EnergyExpressions
using AtomicLevels
using Atoms

makedocs(
    modules = [Atoms],
    sitename = "Atoms",
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Density matrices" => "density_matrices.md"
        ]
    ],
    assets = ["assets/latex.js"]
)

deploydocs(repo = "github.com/JuliaAtoms/Atoms.jl.git")
