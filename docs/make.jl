using Documenter

makedocs(;
    sitename="Norma",
    authors="Norma contributors",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true"),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Running Norma" => "running.md",
        "Testing" => "testing.md",
        "Examples" => "examples.md",
        "Profiling and debugging" => "profiling-and-debugging.md",
        "Troubleshooting" => "troubleshooting.md",
        "Input File Reference" => [
            "Overview" => "reference/index.md",
            "Mesh and output" => "reference/mesh-and-io.md",
            "Model" => "reference/model.md",
            "Materials" => "reference/materials.md",
            "Time integrators" => "reference/time-integrators.md",
            "Solvers" => "reference/solvers.md",
            "Boundary conditions" => "reference/boundary-conditions.md",
            "Initial conditions" => "reference/initial-conditions.md",
            "Function expressions" => "reference/functions.md",
            "Multidomain and Schwarz coupling" => "reference/multidomain.md",
        ],
    ],
    # The reference cross-links with relative .md links; that is intentional and
    # portable, so do not fail the build on the checks that dislike them.
    warnonly=[:cross_references, :missing_docs],
)

deploydocs(;
    repo="github.com/sandialabs/Norma.jl.git",
    devbranch="main",
    push_preview=true,
)
