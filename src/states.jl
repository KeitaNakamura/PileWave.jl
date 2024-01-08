struct ElementState
    # parameters for pile
    E  :: Float64 # young's modulus
    ρ  :: Float64 # density
    A  :: Float64 # diameter
    θ  :: Float64 # perimeter
    # parameters for soil model
    k₁ :: Float64 # stiffness 1
    k₂ :: Float64 # stiffness 2
    τ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield stress factor
    C  :: Float64 # damping
    # variables
    τ̄    :: Float64
    τ̄ₙ   :: Float64
    dτ̄du :: Float64
end

mutable struct ElementStateBottom
    # parameters for pile
    A  :: Float64 # diameter
    # parameters for soil model
    k₁ :: Float64 # stiffness 1
    k₂ :: Float64 # stiffness 2
    σ̄ᵤ :: Float64 # ultimate shaft resistance (friction)
    R  :: Float64 # yield factor
    C  :: Float64 # damping
    # variables
    σ̄    :: Float64
    σ̄ₙ   :: Float64
    dσ̄du :: Float64
    # dof index
    index :: Int
end

function set_elementstate!(est::LazyRow{ElementState}, pile::TOMLPile)
    est.E  = pile.youngs_modulus
    est.ρ  = pile.density
    est.A  = pile.area
    est.θ  = pile.perimeter
    est.k₁ = 0
    est.k₂ = 0
    est.τ̄ᵤ = 0
    est.R  = 0
    est.C  = 0
    est.τ̄    = 0
    est.τ̄ₙ   = 0
    est.dτ̄du = 0
    est
end
function set_elementstate!(est::LazyRow{ElementState}, layer::SoilModel_Shaft)
    est.k₁ = layer.yield_stress / layer.quake_1
    est.k₂ = layer.yield_stress / layer.quake_2
    est.τ̄ᵤ = layer.yield_stress
    est.R  = layer.yield_factor
    est.C  = layer.damping
    est
end
function create_elementstatebottom(pile::TOMLPile, layer::SoilModel_Bottom, btm::Int)
    ElementStateBottom(pile.area_bottom,
                       layer.yield_stress / layer.quake_1,
                       layer.yield_stress / layer.quake_2,
                       layer.yield_stress,
                       layer.yield_factor,
                       layer.damping,
                       0,
                       0,
                       0,
                       btm,)
end

function elementstate_startup!(est::StructArray{ElementState})
    est.τ̄ₙ .= est.τ̄
end
function elementstate_startup!(est::ElementStateBottom)
    est.σ̄ₙ = est.σ̄
end

# vector
function assemble!(
        p::AbstractVector,
        grid::Grid,
        Δu˜::AbstractMatrix,
        v˜::AbstractMatrix,
        est::StructArray{ElementState}
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
        w*θ*(τ̄ + C*v)
    end

    nothing
end
# matrix
function assemble!(
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        grid::Grid,
        Δu˜::AbstractMatrix,
        v˜::AbstractMatrix,
        est::StructArray{ElementState}
    )
    # velocity contribution
    integrate!(C_tan, Sf(), grid) do i, w, dv
        θ = est.θ[i]
        C = est.C[i]
        w*θ*C*dv
    end

    # displacement contribution
    integrate!(K_tan, Sf(), grid) do i, w, du
        θ = est.θ[i]
        dτ̄du = est.dτ̄du[i]
        w*θ*dτ̄du*du
    end

    nothing
end

# vector
function assemble!(
        p::AbstractVector,
        Δu::AbstractVector,
        v::AbstractVector,
        est::ElementStateBottom,
    )
    i = est.index

    # update stress
    dσ̄du, σ̄ = gradient(Δu->staticstress(est, Δu), Δu[i], :all)
    est.σ̄ = σ̄
    est.dσ̄du = dσ̄du

    # assemble
    Aₚ = est.A
    Cₚ = est.C
    p[i] += Aₚ*(σ̄ + Cₚ*v[i])

    nothing
end
# matrix
function assemble!(
        C_tan::AbstractMatrix,
        K_tan::AbstractMatrix,
        Δu::AbstractVector,
        v::AbstractVector,
        est::ElementStateBottom,
    )
    i = est.index

    # update stress
    dσ̄du, σ̄ = gradient(Δu->staticstress(est, Δu), Δu[i], :all)
    est.σ̄ = σ̄

    # assemble
    Aₚ = est.A
    Cₚ = est.C
    dσ̄du = est.dσ̄du
    C_tan[i,i] += Aₚ*Cₚ
    K_tan[i,i] += Aₚ*dσ̄du

    nothing
end

################
# staticstress #
################

function staticstress(k₁::Real, k₂::Real, τy::Real, τₙ::Real, Δu::Real)
    k = ifelse(Δu ≥ 0, k₁, k₂)
    τᵗʳ = τₙ + k * Δu
    fᵗʳ = abs(τᵗʳ) - τy # yield function
    fᵗʳ ≤ 0 && return τᵗʳ
    τᵗʳ - sign(τᵗʳ)*fᵗʳ
end

function staticstress(k₁::Real, k₂::Real, τᵤ::Real, R::Real, τₙ::Real, Δu::Real)
    χ = (τᵤ - R*τᵤ) / 2
    staticstress(k₁, k₂, τᵤ-χ, τₙ-χ, Δu) + χ
end

function staticstress(est::ElementState, Δu::Real)
    staticstress(est.k₁, est.k₂, est.τ̄ᵤ, est.R, est.τ̄ₙ, Δu)
end

function staticstress(est::ElementStateBottom, Δu::Real)
    staticstress(est.k₁, est.k₂, est.σ̄ᵤ, est.R, est.σ̄ₙ, Δu)
end
