module PileWave

using Femto
using TOMLX

using StructArrays
using Interpolations: linear_interpolation

include("input.jl")
include("states.jl")
include("solver.jl")

end # module PileWave
