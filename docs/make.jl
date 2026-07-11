using Documenter

# The experiment figures are produced by the scripts in `scripts/` and written to `plots/`
# (git-ignored). Copy them into the documentation source tree so Documenter includes them in the
# build. Run the experiment scripts first; any figure that is missing simply won't be embedded.
const figdir = joinpath(@__DIR__, "src", "figures")
const plotsdir = joinpath(@__DIR__, "..", "plots")
mkpath(figdir)
if isdir(plotsdir)
    for f in readdir(plotsdir)
        endswith(f, ".png") && cp(joinpath(plotsdir, f), joinpath(figdir, f); force = true)
    end
else
    @warn "no plots/ directory found — run the experiment scripts before building the docs" plotsdir
end

makedocs(;
    sitename = "ReducedPrecision.jl",
    authors  = "Michael Kraus",
    format   = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home"        => "index.md",
        "Methodology" => "methodology.md",
        "Experiments" => [
            "Harmonic Oscillator" => "harmonic_oscillator.md",
            "Pendulum"            => "pendulum.md",
            "Double Pendulum"     => "double_pendulum.md",
            "Toda Lattice"        => "toda_lattice.md",
        ],
        "Findings" => "findings.md",
    ],
)

deploydocs(;
    repo      = "github.com/michakraus/ReducedPrecision.jl.git",
    devbranch = "main",
)
