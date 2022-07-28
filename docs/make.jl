using Documenter
using EnergyExpressions
using AtomicLevels
using Atoms

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(
    modules = [Atoms],
    sitename = "Atoms",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
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
    format = Documenter.HTML(
        mathengine = MathJax2(Dict(:TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(
                :defd => "≝",
                :ket => ["|#1\\rangle",1],
                :bra => ["\\langle#1|",1],
                :braket => ["\\langle#1|#2\\rangle",2],
                :ketbra => ["|#1\\rangle\\!\\langle#2|",2],
                :matrixel => ["\\langle#1|#2|#3\\rangle",3],
                :vec => ["\\mathbf{#1}",1],
                :mat => ["\\mathsf{#1}",1],
                :conj => ["#1^*",1],
                :im => "\\mathrm{i}",
                :operator => ["\\mathfrak{#1}",1],
                :Hamiltonian => "\\operator{H}",
                :hamiltonian => "\\operator{h}",
                :Lagrangian => "\\operator{L}",
                :fock => "\\operator{f}",
                :lagrange => ["\\epsilon_{#1}",1],
                :vary => ["\\delta_{#1}",1],
                :onebody => ["(#1|#2)",2],
                :twobody => ["[#1|#2]",2],
                :twobodydx => ["[#1||#2]",2],
                :direct => ["{\\operator{J}_{#1}}",1],
                :exchange => ["{\\operator{K}_{#1}}",1],
                :diff => ["\\mathrm{d}#1\\,",1],
            ),
        ))),
    ),
    doctest=false,
)

deploydocs(
    repo = "github.com/JuliaAtoms/Atoms.jl.git",
    push_preview = true,
)
