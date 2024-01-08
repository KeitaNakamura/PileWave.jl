const TOMLDict = Dict{String, Any}

##############
# Soil layer #
##############

Base.@kwdef struct SoilModel_Shaft
    quake        :: Float64 = Inf
    quake_1      :: Float64 = quake
    quake_2      :: Float64 = quake_1
    yield_stress :: Float64 = 0.0
    yield_factor :: Float64 = 1.0
    damping      :: Float64 = 0.0
end
Base.@kwdef struct SoilModel_Bottom
    quake        :: Float64 = Inf
    quake_1      :: Float64 = quake
    quake_2      :: Float64 = quake_1
    yield_stress :: Float64 = 0.0
    yield_factor :: Float64 = 0.0
    damping      :: Float64 = 0.0
end

Base.@kwdef struct TOMLSoilLayer
    thickness :: Float64
    shaft     :: SoilModel_Shaft  = SoilModel_Shaft()
    bottom    :: SoilModel_Bottom = SoilModel_Bottom()
end

####################
# Other conditions #
####################

Base.@kwdef struct TOMLInput
    embedded_depth :: Float64
    gravity        :: Float64
    t_stop         :: Float64
    input_load     :: Union{Function, String}
end

Base.@kwdef mutable struct TOMLOutput
    directory      :: String          = ""
    num_data       :: Int             = 200
    paraview       :: Bool            = false
    history_points :: Vector{Float64} = [0.0]
    show_progress  :: Bool            = true
end

Base.@kwdef struct TOMLAdvanced_NewmarkBeta
    beta  :: Float64 = 1/4
    gamma :: Float64 = 1/2
end

Base.@kwdef struct TOMLAdvanced
    shape       :: Shape                    = Line3()
    CFL         :: Float64                  = 0.5
    NewmarkBeta :: TOMLAdvanced_NewmarkBeta = TOMLAdvanced_NewmarkBeta()
end

Base.@kwdef struct TOMLPile
    length         :: Float64
    num_elements   :: Float64
    # material parameters
    youngs_modulus :: Float64
    density        :: Float64
    area           :: Float64
    perimeter      :: Float64
    area_bottom    :: Float64 = area
end

struct TOMLFile
    name       :: String
    Input      :: TOMLInput
    Output     :: TOMLOutput
    Advanced   :: TOMLAdvanced
    Pile       :: Vector{TOMLPile}
    SoilLayer  :: Vector{TOMLSoilLayer}
end

function read_inputfile(file::String)
    dict = TOMLX.parsefile(@__MODULE__, file)
    dict["__name__"] = file
    read_input(dict)
end

function read_input(dict::TOMLDict)
    name = get(dict, "__name__", "fem1d")
    Input      = TOMLX.from_dict(TOMLInput, dict["Input"])
    Output     = TOMLX.from_dict(TOMLOutput, get(dict, "Output", TOMLDict()))
    Advanced   = TOMLX.from_dict(TOMLAdvanced, get(dict, "Advanced", TOMLDict()))
    Pile       = map(pile->TOMLX.from_dict(TOMLPile, pile), dict["Pile"])
    SoilLayer  = map(layer->TOMLX.from_dict(TOMLSoilLayer, layer), dict["SoilLayer"])
    # modify output directory
    Output.directory = Output.directory=="" ? splitext(basename(name))[1]*".tmp" : Output.directory
    TOMLFile(name, Input, Output, Advanced, Pile, SoilLayer)
end
