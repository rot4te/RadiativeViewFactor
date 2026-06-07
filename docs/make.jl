using Documenter
using RadiativeViewFactor

makedocs(
    sitename = "RadiativeViewFactor.jl",
    authors  = "Alex Coxe",
    format   = Documenter.HTML(
        prettyurls       = get(ENV, "CI", nothing) == "true",
        canonical        = "https://rot4te.github.io/RadiativeViewFactor.jl",
        edit_link        = "main",
        assets           = String[],
    ),
    modules  = [RadiativeViewFactor],
    pages    = [
        "Home"           => "index.md",
        "Manual"         => [
            "Getting Started"     => "manual/getting_started.md",
            "Mesh Requirements"   => "manual/mesh_requirements.md",
            "Integration Methods" => "manual/integration_methods.md",
            "Obstruction Detection" => "manual/obstruction.md",
            "GPU Backends"        => "manual/gpu.md",
            "Visualisation"       => "manual/visualisation.md",
            "Performance Guide"   => "manual/performance.md",
        ],
        "Theory"         => "theory.md",
        "API Reference"  => "api.md",
        "References"     => "references.md",
    ],
    checkdocs = :exports,
    warnonly  = true,
)

deploydocs(
    repo   = "github.com/rot4te/RadiativeViewFactor.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
)
