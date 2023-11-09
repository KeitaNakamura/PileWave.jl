using PileWave
using Test

using LinearAlgebra
using TOMLX

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
        @test norm(PileWave.solve(dict; return_u=true)) ≈ ans
    end
end

@testset "FEM 1D" begin
    nlayers = 3
    @testset "Bottom free" begin
        test_fem_multilayers("FEM1D/bottom_free/voigt.toml", nlayers, 35.22443769539338)
        test_fem_multilayers("FEM1D/bottom_free/smith.toml", nlayers, 35.22443769539338)
    end
    @testset "Bottom fixed" begin
        test_fem_multilayers("FEM1D/bottom_fixed/voigt.toml", nlayers, 9.212313943152246)
        test_fem_multilayers("FEM1D/bottom_fixed/smith.toml", nlayers, 9.212313943152246)
    end
    @testset "Bottom friction-dashpot" begin
        test_fem_multilayers("FEM1D/bottom_friction_dashpot/voigt.toml", nlayers, 6.469511528046563)
        test_fem_multilayers("FEM1D/bottom_friction_dashpot/smith.toml", nlayers, 6.388964534822483)
    end
    @testset "Shaft and bottom friction-dashpot" begin
        test_fem_multilayers("FEM1D/shaft_bottom_friction_dashpot/voigt.toml", nlayers, 7.038348438563308)
        test_fem_multilayers("FEM1D/shaft_bottom_friction_dashpot/smith.toml", nlayers, 6.818677234812194)
    end
    @testset "Complex version" begin
        @test norm(PileWave.solve("FEM1D/complex/voigt.toml"; return_u=true)) ≈ 0.17796429382389684
        @test norm(PileWave.solve("FEM1D/complex/smith.toml"; return_u=true)) ≈ 0.12515693582582205
    end
end
