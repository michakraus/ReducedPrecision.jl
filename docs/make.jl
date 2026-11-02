using Documenter

# The experiment figures are produced by the scripts in `scripts/`, which write them to `plots/`
# (git-ignored). In CI one job per script runs it and the `documenter` job collects all of them
# into `plots/`; locally, run `julia --project=. scripts/run_all.jl` first. They are copied into
# the documentation source tree, so that Documenter includes them in the build.
const srcdir = joinpath(@__DIR__, "src")
const figdir = joinpath(srcdir, "figures")
const plotsdir = joinpath(@__DIR__, "..", "plots")

isdir(plotsdir) ||
    error("no plots/ directory at $(normpath(plotsdir)): run the experiment scripts in " *
          "scripts/ before building the docs")

mkpath(figdir)
for f in readdir(plotsdir)
    endswith(f, ".png") && cp(joinpath(plotsdir, f), joinpath(figdir, f); force = true)
end

# Fail before `makedocs` if a page links a figure that no script produced, and name the
# figures and the pages. The filename prefix is the script's problem name, and the larger of a
# problem's two time steps comes from its `_longtime` variant. The check reads `plots/`, not
# `figures/`, which keeps the copies of earlier builds.
missing_figures = String[]
for page in filter(endswith(".md"), readdir(srcdir))
    for m in eachmatch(r"\]\((figures/[^)\s]+)\)", read(joinpath(srcdir, page), String))
        isfile(joinpath(plotsdir, basename(m[1]))) ||
            push!(missing_figures, "$page: $(m[1])")
    end
end
isempty(missing_figures) ||
    error("$(length(missing_figures)) figure(s) linked from docs/src are missing from " *
          "plots/:\n  " * join(missing_figures, "\n  "))

makedocs(;
    sitename = "ReducedPrecision.jl",
    authors = "Michael Kraus",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Methodology" => "methodology.md",
        "Experiments" => [
            "Harmonic Oscillator" => "harmonic_oscillator.md",
            "Pendulum" => "pendulum.md",
            "Double Pendulum" => "double_pendulum.md",
            "Toda Lattice" => "toda_lattice.md",
            "Lotka–Volterra 2D" => "lotka_volterra_2d.md",
            "Lotka–Volterra 4D" => "lotka_volterra_4d.md"
        ],
        "Findings" => "findings.md"
    ]
)

# Publish a complete site only. `DEPLOY_DOCS` is how the workflow withholds a deployment when
# some script job failed. Defaults to deploying, so that local builds and manual runs are
# unaffected.
if get(ENV, "DEPLOY_DOCS", "true") == "true"
    deploydocs(;
        repo = "github.com/michakraus/ReducedPrecision.jl.git",
        devbranch = "main"
    )
else
    @warn "DEPLOY_DOCS is not true: some figure job failed, so the site is not deployed"
end
