abstract type ElementState end
abstract type ElementStateBottom end

################
# staticstress #
################

function staticstress(q₁::Real, q₂::Real, τᵤ::Real, R::Real, τₙ::Real, Δu::Real)
    (k, τy) = Δu ≥ 0 ? (τᵤ/q₁, τᵤ) : (τᵤ/q₂, R*τᵤ)
    τᵗʳ = τₙ + k * Δu
    fᵗʳ = abs(τᵗʳ) - τy # yield function
    if fᵗʳ ≤ 0
        return τᵗʳ
    else
        return τᵗʳ - sign(τᵗʳ)*fᵗʳ
    end
end

function staticstress(est::ElementState, Δu::Real)
    staticstress(est.q₁, est.q₂, est.τ̄ᵤ, est.R, est.τ̄ₙ, Δu)
end

function staticstress(est::ElementStateBottom, Δu::Real)
    staticstress(est.q₁, est.q₂, est.σ̄ᵤ, est.R, est.σ̄ₙ, Δu)
end

###############
# Voigt model #
###############

struct VoigtElementState <: ElementState
    # parameters for pile
    E  :: Float64 # young's modulus
    ρ  :: Float64 # density
    A  :: Float64 # diameter
    θ  :: Float64 # perimeter
    # parameters for soil model
    q₁ :: Float64 # quake 1
    q₂ :: Float64 # quake 2
    τ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield stress factor
    C  :: Float64 # damping
    # variables
    τ̄    :: Float64
    τ̄ₙ   :: Float64
    dτ̄du :: Float64
end

mutable struct VoigtElementStateBottom <: ElementStateBottom
    # parameters for pile
    A  :: Float64 # diameter
    # parameters for soil model
    q₁ :: Float64 # quake 1
    q₂ :: Float64 # quake 2
    σ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield factor
    C  :: Float64 # damping
    # variables
    σ̄    :: Float64
    σ̄ₙ   :: Float64
    # dof index
    index :: Int
end

get_elementstate_type(::Type{VoigtModel}) = VoigtElementState
function set_elementstate!(est::LazyRow{VoigtElementState}, pile::TOMLPile)
    est.E  = pile.youngs_modulus
    est.ρ  = pile.density
    est.A  = pile.area
    est.θ  = pile.perimeter
    est.q₁ = Inf
    est.q₂ = Inf
    est.τ̄ᵤ = 0
    est.R  = 0
    est.C  = 0
    est.τ̄    = 0
    est.τ̄ₙ   = 0
    est.dτ̄du = 0
    est
end
function set_elementstate!(est::LazyRow{VoigtElementState}, layer::VoigtModel)
    est.q₁ = layer.quake_1
    est.q₂ = layer.quake_2
    est.τ̄ᵤ = layer.yield_stress
    est.R  = layer.yield_factor
    est.C  = layer.damping
    est
end
function create_elementstatebottom(::Type{VoigtModel}, pilebottom::TOMLPileBottom, layer::VoigtModel, btm::Int)
    VoigtElementStateBottom(pilebottom.area,
                            layer.quake_bottom_1,
                            layer.quake_bottom_2,
                            layer.yield_stress_bottom,
                            layer.yield_factor_bottom,
                            layer.damping_bottom,
                            0,
                            0,
                            btm,)
end

function elementstate_startup!(est::StructArray{VoigtElementState})
    est.τ̄ₙ .= est.τ̄
end
function elementstate_startup!(est::VoigtElementStateBottom)
    est.σ̄ₙ = est.σ̄
end

function assemble!(
        p::AbstractVector,
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        grid::Grid,
        Δu˜::AbstractMatrix,
        v˜::AbstractMatrix,
        est::StructArray{VoigtElementState}
    )
    # update stresses
    for i in eachindex(est)
        Δu = Δu˜[i]
        dτ̄du, τ̄ = gradient(Δu->staticstress(est[i], Δu), Δu, :all)
        est.τ̄[i] = τ̄
        est.dτ̄du[i] = dτ̄du
    end

    # force
    integrate!(p, Sf(), grid) do i, w
        θ = est.θ[i]
        C = est.C[i]
        τ̄ = est.τ̄[i]
        v = v˜[i]
        w*(τ̄*θ + C*v)
    end

    # velocity contribution
    integrate!(C_tan, Sf(), grid) do i, w, dv
        C = est.C[i]
        w*C*dv
    end

    # displacement contribution
    integrate!(K_tan, Sf(), grid) do i, w, du
        θ = est.θ[i]
        dτ̄du = est.dτ̄du[i]
        w*dτ̄du*θ*du
    end

    nothing
end

function assemble!(
        p::AbstractVector,
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        Δu::AbstractVector,
        v::AbstractVector,
        est::VoigtElementStateBottom,
    )
    i = est.index

    # update stress
    dσ̄du, σ̄ = gradient(Δu->staticstress(est, Δu), Δu[i], :all)
    est.σ̄ = σ̄

    # assemble
    Aₚ = est.A
    Cₚ = est.C
    p[i] += Aₚ*(σ̄ + Cₚ*v[i])
    C_tan[i,i] += Aₚ*Cₚ
    K_tan[i,i] += Aₚ*dσ̄du

    nothing
end

###############
# Smith model #
###############

struct SmithElementState <: ElementState
    # parameters for pile
    E  :: Float64 # young's modulus
    ρ  :: Float64 # density
    A  :: Float64 # diameter
    θ  :: Float64 # perimeter
    # parameters for soil model
    q₁ :: Float64 # quake 1
    q₂ :: Float64 # quake 2
    τ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield stress factor
    J  :: Float64 # smith damping
    # variables
    τ̄    :: Float64
    τ̄ₙ   :: Float64
    dτ̄du :: Float64
end

mutable struct SmithElementStateBottom <: ElementStateBottom
    # parameters for pile
    A  :: Float64 # diameter
    # parameters for soil model
    q₁ :: Float64 # quake 1
    q₂ :: Float64 # quake 2
    σ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield factor
    J  :: Float64 # smith damping
    # variables
    σ̄    :: Float64
    σ̄ₙ   :: Float64
    # dof index
    index :: Int
end

get_elementstate_type(::Type{SmithModel}) = SmithElementState
function set_elementstate!(est::LazyRow{SmithElementState}, pile::TOMLPile)
    est.E  = pile.youngs_modulus
    est.ρ  = pile.density
    est.A  = pile.area
    est.θ  = pile.perimeter
    est.q₁ = Inf
    est.q₂ = Inf
    est.τ̄ᵤ = 0
    est.R  = 0
    est.J  = 0
    est.τ̄    = 0
    est.τ̄ₙ   = 0
    est.dτ̄du = 0
    est
end
function set_elementstate!(est::LazyRow{SmithElementState}, layer::SmithModel)
    est.q₁ = layer.quake_1
    est.q₂ = layer.quake_2
    est.τ̄ᵤ = layer.yield_stress
    est.R  = layer.yield_factor
    est.J  = layer.damping
    est
end
function create_elementstatebottom(::Type{SmithModel}, pilebottom::TOMLPileBottom, layer::SmithModel, btm::Int)
    SmithElementStateBottom(pilebottom.area,
                            layer.quake_bottom_1,
                            layer.quake_bottom_2,
                            layer.yield_stress_bottom,
                            layer.yield_factor_bottom,
                            layer.damping_bottom,
                            0,
                            0,
                            btm,)
end

function elementstate_startup!(est::StructArray{SmithElementState})
    est.τ̄ₙ .= est.τ̄
end
function elementstate_startup!(est::SmithElementStateBottom)
    est.σ̄ₙ = est.σ̄
end

function assemble!(
        p::AbstractVector,
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        grid::Grid,
        Δu˜::AbstractMatrix,
        v˜::AbstractMatrix,
        est::StructArray{SmithElementState}
    )
    # update stresses
    for i in eachindex(est)
        Δu = Δu˜[i]
        dτ̄du, τ̄ = gradient(Δu->staticstress(est[i], Δu), Δu, :all)
        est.τ̄[i] = τ̄
        est.dτ̄du[i] = dτ̄du
    end

    # force
    integrate!(p, Sf(), grid) do i, w
        θ = est.θ[i]
        J = est.J[i]
        τ̄ = est.τ̄[i]
        v = v˜[i]
        w*τ̄*θ*(1 + J*v)
    end

    # velocity contribution
    integrate!(C_tan, Sf(), grid) do i, w, dv
        θ = est.θ[i]
        J = est.J[i]
        τ̄ = est.τ̄[i]
        w*τ̄*θ*J*dv
    end

    # displacement contribution
    integrate!(K_tan, Sf(), grid) do i, w, du
        θ = est.θ[i]
        J = est.J[i]
        dτ̄du = est.dτ̄du[i]
        v = v˜[i]
        w*dτ̄du*θ*(1 + J*v)*du
    end

    nothing
end

function assemble!(
        p::AbstractVector,
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        Δu::AbstractVector,
        v::AbstractVector,
        est::SmithElementStateBottom,
    )
    i = est.index

    # update stress
    dσ̄du, σ̄ = gradient(Δu->staticstress(est, Δu), Δu[i], :all)
    est.σ̄ = σ̄

    # assemble
    Aₚ = est.A
    Jₚ = est.J
    p[i] += σ̄*Aₚ*(1 + Jₚ*v[i])
    C_tan[i,i] += σ̄*Aₚ*Jₚ
    K_tan[i,i] += dσ̄du*Aₚ*(1 + Jₚ*v[i])

    nothing
end
