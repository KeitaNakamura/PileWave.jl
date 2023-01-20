module FEM

using Femto
using TOMLX

using StructArrays
using Interpolations: linear_interpolation

include("input.jl")
include("states.jl")
include("core.jl")

end # module FEM
