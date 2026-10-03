using Documenter
using DocumenterVitepress
using BiophysicalBehaviour
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

# Drawings of bodies, from the documentation of BiophysicalGeometry.jl, loaded with `using Main.GeometryFigures`
include("geometry_figures.jl")

# The microclimates of "Get started", loaded in the manual with `using Main.ExampleEnvironments`
include("example_environments.jl")

makedocs(
    modules = [BiophysicalBehaviour],
    sitename = "BiophysicalBehaviour.jl",
    authors = "Michael Kearney, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/BiophysicalBehaviour.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/BiophysicalBehaviour.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
