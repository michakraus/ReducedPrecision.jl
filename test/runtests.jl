using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "BFloat16 compatibility shims" include("bfloat16_compat.jl")
    @safetestset "Method registry" include("methods.jl")
    @safetestset "Tableau-driven initial guess" include("initial_guess.jl")
    @safetestset "Study" include("study.jl")
    @safetestset "Diagnostics" include("diagnostics.jl")
    @safetestset "Plotting" include("plotting.jl")
end
