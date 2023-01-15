###############
# Soil models #
###############

abstract type SoilModel end

Base.@kwdef struct VoigtModel <: SoilModel
    thickness           :: Float64
    # shaft
    quake_1             :: Float64
    quake_2             :: Float64 = quake_1
    yield_stress        :: Float64
    yield_factor        :: Float64 = 1.0
    damping             :: Float64
    # bottom
    quake_bottom_1      :: Float64 = Inf
    quake_bottom_2      :: Float64 = quake_bottom_1
    yield_stress_bottom :: Float64 = 0.0
    yield_factor_bottom :: Float64 = 0.0
    damping_bottom      :: Float64 = 0.0
end

Base.@kwdef struct SmithModel <: SoilModel
    thickness           :: Float64
    # shaft
    quake_1             :: Float64
    quake_2             :: Float64 = quake_1
    yield_stress        :: Float64
    yield_factor        :: Float64 = 1.0
    damping             :: Float64
    # bottom
    quake_bottom_1      :: Float64 = Inf
    quake_bottom_2      :: Float64 = quake_bottom_1
    yield_stress_bottom :: Float64 = 0.0
    yield_factor_bottom :: Float64 = 0.0
    damping_bottom      :: Float64 = 0.0
end

####################
# Other conditions #
####################

Base.@kwdef struct TOMLInput{Model}
    soilmodel      :: Type{Model}
    embedded_depth :: Float64
    gravity        :: Float64 = 0.0
    t_stop         :: Float64
    loadinput      :: Union{Function, String}
end

Base.@kwdef struct TOMLOutput
    directory      :: String      = "Output"
    num_data       :: Int         = 100
    history_points :: Vector{Int} = Int[]
end

Base.@kwdef struct TOMLAdvanced{S <: Shape}
    shape :: S       = Line3()
    CFL   :: Float64 = 1.0
end

Base.@kwdef struct TOMLPile
    length         :: Float64
    num_elements   :: Float64
    # material parameters
    youngs_modulus :: Float64
    density        :: Float64
    area           :: Float64
    perimeter      :: Float64
end

Base.@kwdef struct TOMLPileBottom
    area :: Float64
end

struct TOMLFile{M <: SoilModel, S <: Femto.Line}
    name       :: String
    Input      :: TOMLInput{M}
    Output     :: TOMLOutput
    Advanced   :: TOMLAdvanced{S}
    Pile       :: Vector{TOMLPile}
    PileBottom :: TOMLPileBottom
    SoilLayer  :: Vector{M}
end

function read_inputfile(file::String)
    dict = TOMLX.parsefile(@__MODULE__, file)
    Input      = TOMLX.from_dict(TOMLInput, dict["Input"])
    Output     = TOMLX.from_dict(TOMLOutput, get(dict, "Output", Dict{String, Any}()))
    Advanced   = TOMLX.from_dict(TOMLAdvanced, get(dict, "Advanced", Dict{String, Any}()))
    Pile       = map(pile->TOMLX.from_dict(TOMLPile, pile), dict["Pile"])
    PileBottom = TOMLX.from_dict(TOMLPileBottom, dict["PileBottom"])
    SoilLayer  = map(layer->TOMLX.from_dict(Input.soilmodel, layer), dict["SoilLayer"])
    TOMLFile(file, Input, Output, Advanced, Pile, PileBottom, SoilLayer)
end
