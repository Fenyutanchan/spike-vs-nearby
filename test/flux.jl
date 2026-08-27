# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

@testset "isotropic density-to-flux relation" begin
    density = EU(8.0, 2)

    @test unit_isapprox(
        isotropic_differential_flux(density),
        density / (4 * pi),
    )
    @test unit_isapprox(
        isotropic_differential_flux(density; beta=0.5),
        density / (8 * pi),
    )
end
