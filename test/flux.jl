# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using NaturalUnits: GeV, MeV

@testset "Experimental differential flux units" begin
    @test unit_isapprox(differential_flux_unit, inv(GeV(1.0) * NU.m^2 * NU.s))
    @test flux_value(2.5 * differential_flux_unit) ≈ 2.5
    @test flux_value(convert(MeV, 2.5 * differential_flux_unit)) ≈ 2.5
    # Flux conversion must not silently discard a physical dimension.
    @test_throws DimensionMismatch flux_value(GeV(1.0))
end

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
