module PileWave

using Femto
using TOMLX
using CSV

using StructArrays
using Interpolations: linear_interpolation

include("input.jl")
include("states.jl")
include("solver.jl")

function julia_main()::Cint
    if isempty(ARGS)
        inputtoml = "input.toml"
    else
        inputtoml = ARGS[1]
    end
    try
        solve(inputtoml)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

end # module PileWave
