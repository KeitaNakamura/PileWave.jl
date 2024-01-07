using LinearAlgebra
using TOML
const TOMLDict = Dict{String, Any}
# function test_fem_multilayers(file::String, nlayers::Int, ans::Real)
    # for n in 1:nlayers
        # dict = TOMLX.parsefile(PileWave, file)
        # sl1 = dict["SoilLayer"][begin]
        # t = sl1["thickness"]
        # dict["SoilLayer"] = map(1:n) do i
            # sl = copy(sl1)
            # sl["thickness"] = t / n
            # sl
        # end
        # dict["__name__"] = file
        # @test norm(PileWave.solve(dict; return_solution=true).u[:,end]) ≈ ans
    # end
# end

# @testset "Examples" begin
    # nlayers = 3
    # @testset "Bottom free" begin
        # test_fem_multilayers("FEM1D/bottom_free/voigt.toml", nlayers, 35.22443769539338)
        # test_fem_multilayers("FEM1D/bottom_free/smith.toml", nlayers, 35.22443769539338)
    # end
    # @testset "Bottom fixed" begin
        # test_fem_multilayers("FEM1D/bottom_fixed/voigt.toml", nlayers, 9.212313943152246)
        # test_fem_multilayers("FEM1D/bottom_fixed/smith.toml", nlayers, 9.212313943152246)
    # end
    # @testset "Bottom friction-dashpot" begin
        # test_fem_multilayers("FEM1D/bottom_friction_dashpot/voigt.toml", nlayers, 6.469511528046563)
        # test_fem_multilayers("FEM1D/bottom_friction_dashpot/smith.toml", nlayers, 6.388964534822483)
    # end
    # @testset "Shaft and bottom friction-dashpot" begin
        # test_fem_multilayers("FEM1D/shaft_bottom_friction_dashpot/voigt.toml", nlayers, 7.038348438563308)
        # test_fem_multilayers("FEM1D/shaft_bottom_friction_dashpot/smith.toml", nlayers, 6.818677234812194)
    # end
    # @testset "Complex version" begin
        # @test norm(PileWave.solve("FEM1D/complex/voigt.toml"; return_solution=true).u[:,end]) ≈ 0.12592941314005968
        # @test norm(PileWave.solve("FEM1D/complex/smith.toml"; return_solution=true).u[:,end]) ≈ 0.04020891947776435
    # end
# end

@testset "Analytical solutions" begin
    for SoilModel in (PileWave.VoigtModel, PileWave.SmithModel)
        R = 200e-3
        L = 10.0
        E = 3e10
        A = π*R^2
        θ = 2π*R
        ρ = 1.2e3
        F = 20e6 * A
        default = merge(
            TOMLDict(
                "Input" => TOMLDict(
                    "soilmodel"      => SoilModel,
                    "embedded_depth" => L,
                    "gravity"        => 0.0,
                    "t_stop"         => 4e-3,
                    "input_load"     => t -> ifelse(t<1e-3, F*sin(π*(t/1e-3)), 0)
                ),
                # "Output" => TOMLDict(
                    # "num_data" => 500,
                # ),
                "Advanced" => TOMLDict(
                    "CFL" => 0.1,
                ),
            ),
            TOML.parse("""
            [[Pile]]
            length         = $L
            num_elements   = 100
            youngs_modulus = $E
            density        = $ρ
            area           = $A
            perimeter      = $θ
            """)
        )
        @testset "Simple sine wave with no friction on shaft" begin
            fcurve(z) = z < 5 ? F*sin(π*(z/5)) : 0.0
            @testset "Bottom fixed" begin
                dict = merge(default, TOML.parse("""
                [[SoilLayer]]
                thickness = 20.0
                bottom = {yield_stress = inf}
                """))
                sol = PileWave.solve(dict, return_solution=true)
                @test sol.f[:,end] ≈ map(fcurve, sol.z) rtol=0.01
            end
            @testset "Bottom free" begin
                dict = merge(default, TOML.parse("""
                [[SoilLayer]]
                thickness = 20.0
                """))
                sol = PileWave.solve(dict, return_solution=true)
                @test sol.f[:,end] ≈ -map(fcurve, sol.z) rtol=0.01
            end
        end
        @testset "Vibration motion" begin
            default["Advanced"]["CFL"] = 10.0
            @testset "Spring on shaft" begin
                default["Input"]["t_stop"] = 100e-3
                k = 2e6
                σ_y = 100e6 # should be large enough
                dict = merge(default, TOML.parse("""
                [[SoilLayer]]
                thickness = 20.0
                shaft = {quake = $(σ_y/k), yield_stress = $σ_y}
                """))
                sol = PileWave.solve(dict, return_solution=true)
                index = argmin(i->abs(sol.z[i]-5), eachindex(sol.z))
                a = maximum(sol.u[index,:])
                m = ρ*A*L
                K = k*θ*L
                ω = √(K/m)
                @test sol.u[index,:] ≈ map(t->a*sin(ω*t), sol.t) rtol=0.12
            end
            if SoilModel isa PileWave.VoigtModel
                @testset "Spring and dashpot on shaft" begin
                    default["Input"]["t_stop"] = 200e-3
                    k = 2e6
                    σ_y = 1000e6 # should be large enough
                    c = 5e3
                    dict = merge(default, TOML.parse("""
                    [[SoilLayer]]
                    thickness = 20.0
                    shaft = {quake = $(σ_y/k), yield_stress = $σ_y, damping = $c}
                    """))
                    sol = PileWave.solve(dict, return_solution=true)
                    index = argmin(i->abs(sol.z[i]-5), eachindex(sol.z))
                    m = ρ*A*L
                    K = k*θ*L
                    C = c*θ*L
                    ω = 1/2m*√(4m*K-C^2)
                    Γ = C/2m
                    T = 2π/ω
                    a = maximum(sol.u[index,:]) * exp(Γ*T/4)
                    @test sol.u[index,:] ≈ map(t->a*exp(-Γ*t)*sin(ω*t), sol.t) rtol=0.12
                end
            end
        end
    end
end

@testset "Nearly static loading" begin
    for SoilModel in (PileWave.VoigtModel, PileWave.SmithModel)
        F = 500.0
        t1 = 10.0 # linearly incrase the load by `t1`
        L = 10.0
        E = 10e3
        A = 1.0
        θ = 1.0
        default = merge(
            TOMLDict(
                "Input" => TOMLDict(
                    "soilmodel" => SoilModel,
                    "embedded_depth" => 10.0,
                    "gravity" => 0.0,
                    "t_stop" => 500.0,
                    "input_load" => t -> ifelse(t < t1, (F/t1)*t, F)
                ),
                "Advanced" => TOMLDict(
                    "shape" => PileWave.Line2(),
                    "CFL" => 10.0,
                ),
            ),
            TOML.parse("""
            [[Pile]]
            length         = $L
            num_elements   = 100
            youngs_modulus = $E
            density        = 1.0
            area           = $A
            perimeter      = $θ
            """)
        )
        @testset "No friction on shaft" begin
            # bottom fixed
            dict = merge(default, TOML.parse("""
            [[SoilLayer]]
            thickness = 20.0
            bottom = {yield_stress = inf}
            """))
            sol = PileWave.solve(dict, return_solution=true)
            u_top = F/(E*A) * L
            @test sol.u[:,end] ≈ LinRange(u_top, 0, 101) rtol=0.01
            @test sol.f[:,end] ≈ LinRange(F, F, 101) rtol=0.01

            # spring at the bottom
            q = 1.0
            σ_y = 2e3 # should be large enough
            dict = merge(default, TOML.parse("""
            [[SoilLayer]]
            thickness = 20.0
            bottom = {quake = $q, yield_stress = $σ_y}
            """))
            sol = PileWave.solve(dict, return_solution=true)
            k = σ_y / q
            u_bottom = F / k
            @test sol.u[:,end] ≈ LinRange(u_top+u_bottom, u_bottom, 101) rtol=0.01
            @test sol.f[:,end] ≈ LinRange(F, F, 101) rtol=0.01
        end
        @testset "Friction on shaft" begin
            σ_y = 10.0

            # bottom fixed
            dict = merge(default, TOML.parse("""
            [[SoilLayer]]
            thickness = 20.0
            shaft = {quake = 1e-3, yield_stress = $σ_y} # immediately slip
            bottom = {yield_stress = inf}
            """))
            f_resist = σ_y * L # we can assume this because of the small quake
            f_top = F
            ϵ_top = f_top / E
            f_bottom = f_top - f_resist
            ϵ_bottom = f_bottom / E
            u_top = (ϵ_top + ϵ_bottom) * L / 2
            sol = PileWave.solve(dict, return_solution=true)
            @test sol.u[:,end] ≈ LinRange(u_top, 0, 101) rtol=0.05
            @test sol.f[:,end] ≈ LinRange(f_top, f_bottom, 101) rtol=0.05

            # spring at the bottom
            q_bottom = 1.0
            σ_y_bottom = 2e3 # should be large enough
            dict = merge(default, TOML.parse("""
            [[SoilLayer]]
            thickness = 20.0
            shaft = {quake = 1e-3, yield_stress = $σ_y} # immediately slip
            bottom = {quake = $q_bottom, yield_stress = $σ_y_bottom}
            """))
            k_bottom = σ_y_bottom / q_bottom
            u_bottom = f_bottom / k_bottom
            sol = PileWave.solve(dict, return_solution=true)
            @test sol.u[:,end] ≈ LinRange(u_top+u_bottom, u_bottom, 101) rtol=0.05
            @test sol.f[:,end] ≈ LinRange(f_top, f_bottom, 101) rtol=0.05
        end
    end
end
