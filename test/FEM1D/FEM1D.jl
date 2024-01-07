function test_fem_multilayers(file::String, nlayers::Int, ans::Real)
    for n in 1:nlayers
        dict = TOMLX.parsefile(PileWave, file)
        sl1 = dict["SoilLayer"][begin]
        t = sl1["thickness"]
        dict["SoilLayer"] = map(1:n) do i
            sl = copy(sl1)
            sl["thickness"] = t / n
            sl
        end
        dict["__name__"] = file
        @test norm(PileWave.solve(dict; return_solution=true).u[:,end]) ≈ ans
    end
end

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

@testset "Compare with analytical solutions" begin
    for SoilModel in (PileWave.VoigtModel, PileWave.SmithModel)
        @testset "No friction on shaft" begin
            R = 400e-3
            L = 10.0
            E = 3e10
            A = π*R^2
            θ = 2π*R
            ρ = 1.2e3
            default = Dict{String, Any}(
                "Input" => Dict{String, Any}(
                    "soilmodel" => SoilModel,
                    "embedded_depth" => L,
                    "gravity" => 0.0,
                    "t_stop" => 4e-3,
                    "input_load" => t -> ifelse(t<1e-3, 20e6*sin(π*(t/1e-3)), 0)
                ),
                "Advanced" => Dict{String, Any}(
                    "CFL" => 0.1,
                ),
                "Pile" => [
                    Dict{String, Any}(
                        "length" => L,
                        "num_elements" => 100,
                        "youngs_modulus" => E,
                        "density" => ρ,
                        "area" => A,
                        "perimeter" => θ,
                    )
                ],
            )
            @testset "Bottom fixed" begin
                dict = merge(
                    default,
                    Dict{String, Any}(
                        "SoilLayer" => [
                            Dict{String, Any}(
                                "thickness" => 20.0,
                                "bottom" => Dict{String, Any}(
                                    "yield_stress" => Inf,
                                ),
                            )
                        ],
                    )
                )
                sol = PileWave.solve(dict, return_solution=true)
                @test sol.f[:,end] ≈ map(z -> z<5 ? 20e6*sin(π*(z/5)) : 0, sol.z) rtol=1e-2
            end
            @testset "Bottom free" begin
                dict = merge(
                    default,
                    Dict{String, Any}(
                        "SoilLayer" => [
                            Dict{String, Any}(
                                "thickness" => 20.0,
                            )
                        ],
                    )
                )
                sol = PileWave.solve(dict, return_solution=true)
                @test sol.f[:,end] ≈ map(z -> z<5 ? -20e6*sin(π*(z/5)) : 0, sol.z) rtol=1e-2
            end
        end
    end
end

@testset "Check with static behavior" begin
    for SoilModel in (PileWave.VoigtModel, PileWave.SmithModel)
        F = 500.0
        t1 = 10.0 # linearly incrase the load by `t1`
        L = 10.0
        E = 10e3
        A = 1.0
        θ = 1.0

        default = Dict{String, Any}(
            "Input" => Dict{String, Any}(
                "soilmodel" => SoilModel,
                "embedded_depth" => 10.0,
                "gravity" => 0.0,
                "t_stop" => 500.0,
                "input_load" => t -> ifelse(t < t1, (F/t1)*t, F)
            ),
            "Advanced" => Dict{String, Any}(
                "shape" => PileWave.Line2(),
                "CFL" => 10.0,
            ),
            "Pile" => [
                Dict{String, Any}(
                    "length" => L,
                    "num_elements" => 100,
                    "youngs_modulus" => E,
                    "density" => 1.0,
                    "area" => A,
                    "perimeter" => θ,
                )
            ],
        )

        @testset "No friction on shaft" begin
            # bottom fixed
            dict = merge(
                default,
                Dict{String, Any}(
                    "SoilLayer" => [
                        Dict{String, Any}(
                            "thickness" => 20.0,
                            "bottom" => Dict{String, Any}(
                                "yield_stress" => Inf,
                            ),
                        )
                    ],
                )
            )
            u = PileWave.solve(dict, return_solution=true).u[:,end]
            u_top = F/(E*A) * L
            @test u ≈ LinRange(u_top, 0, 101) rtol=1e-2

            # apply bottom condition
            q = 1.0
            σ_y = 2e3
            dict = merge(
                default,
                Dict{String, Any}(
                    "SoilLayer" => [
                        Dict{String, Any}(
                            "thickness" => 20.0,
                            "bottom" => Dict{String, Any}(
                                "quake" => q,
                                "yield_stress" => σ_y,
                            ),
                        )
                    ],
                )
            )
            u = PileWave.solve(dict, return_solution=true).u[:,end]
            k = σ_y / q
            u_bottom = F / k
            @test u ≈ LinRange(u_top+u_bottom, u_bottom, 101) rtol=1e-2
        end

        @testset "Friction on shaft" begin
            σ_y = 10.0

            # bottom fixed
            dict = merge(
                default,
                Dict{String, Any}(
                    "SoilLayer" => [
                        Dict{String, Any}(
                            "thickness" => 20.0,
                            "shaft" => Dict{String, Any}(
                                "quake" => 1e-3, # immediately slip
                                "yield_stress" => σ_y,
                            ),
                            "bottom" => Dict{String, Any}(
                                "yield_stress" => Inf,
                            ),
                        )
                    ],
                )
            )
            f_resist = σ_y * L # we can assume this because of the small quake
            f_top = F
            ϵ_top = f_top / E
            f_bottom = f_top - f_resist
            ϵ_bottom = f_bottom / E
            u_top = (ϵ_top + ϵ_bottom) * L / 2
            u = PileWave.solve(dict, return_solution=true).u[:,end]
            @test u ≈ LinRange(u_top, 0, 101) rtol=0.05

            # apply bottom condition
            q_bottom = 1.0
            σ_y_bottom = 2e3
            dict = merge(
                default,
                Dict{String, Any}(
                    "SoilLayer" => [
                        Dict{String, Any}(
                            "thickness" => 20.0,
                            "shaft" => Dict{String, Any}(
                                "quake" => 1e-3, # immediately slip
                                "yield_stress" => σ_y,
                            ),
                            "bottom" => Dict{String, Any}(
                                "quake" => q_bottom,
                                "yield_stress" => σ_y_bottom,
                            ),
                        )
                    ],
                )
            )
            k_bottom = σ_y_bottom / q_bottom
            u_bottom = f_bottom / k_bottom
            u = PileWave.solve(dict, return_solution=true).u[:,end]
            @test u ≈ LinRange(u_top+u_bottom, u_bottom, 101) rtol=0.05
        end
    end
end
