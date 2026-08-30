# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

@testset "ATNF pulsar catalogue" begin
    @test Geminga_ATNF.name == "J0633+1746"
    @test Geminga_ATNF["RAJ"].value == "06:33:54.1530"
    @test Geminga_ATNF["RAJ"].uncertainty == "28"
    @test Geminga_ATNF["RAJ"].unit == "hms"
    @test Geminga_ATNF["RAJ"].reference == "clm+98"
    @test Geminga_ATNF["DIST_A"].value == "0.19"
    @test Geminga_ATNF["DIST_A"].uncertainty == "7"
    @test Geminga_ATNF["DIST_A"].unit == "kpc"
    @test Geminga_ATNF["DIST_A"].reference == "wan11"

    @test unit_isapprox(
        pulsar_characteristic_age(Geminga_ATNF),
        342 * kyr;
        rtol=2e-3,
    )
    @test unit_isapprox(
        pulsar_spin_down_luminosity(Geminga_ATNF),
        3.25e34 * erg / NU.s;
        rtol=2e-3,
    )

    pulsar = read_ATNF_pulsar("B0011+47")
    @test pulsar.name == "J0014+4746"
    @test pulsar["PSRB"].value == "B0011+47"
end

@testset "AMS-02 electron and positron fluxes" begin
    electrons = read_AMS02_electron_flux()
    positrons = read_AMS02_positron_flux()
    flux_unit = inv(one(EU) * NU.m^2 * NU.s)
    combined = combine_AMS02_electron_positron_flux(
        electrons,
        positrons;
        systematic_correlation=0,
    )
    fully_correlated = combine_AMS02_electron_positron_flux(
        electrons,
        positrons;
        systematic_correlation=1,
    )

    @test length(electrons) == 75
    @test length(positrons) == length(combined) == 74
    @test unit_isapprox(combined[1].energy, EU(0.57))
    @test unit_isapprox(
        combined[1].flux,
        19.979 * flux_unit,
    )
    @test unit_isapprox(
        combined[1].statistical_error,
        hypot(0.021e1, 0.109e0) * flux_unit,
    )
    @test unit_isapprox(
        combined[1].systematic_error,
        hypot(0.078e1, 0.141e0) * flux_unit,
    )
    @test unit_isapprox(
        fully_correlated[1].systematic_error,
        (0.078e1 + 0.141e0) * flux_unit,
    )
    @test unit_isapprox(combined[end].energy_max, EU(1000))
    @test unit_isapprox(
        combined[end].flux,
        (1.771e-7 + 1.927e-8) * flux_unit,
    )
end

@testset "DAMPE electron and positron flux" begin
    measurements = read_DAMPE_electron_positron_flux()

    @test length(measurements) == 38
    @test unit_isapprox(measurements[1].energy, EU(25.7))
    @test unit_isapprox(
        measurements[1].flux,
        1.16e-2 * inv(one(EU) * NU.m^2 * NU.s),
    )
    @test unit_isapprox(measurements[end].energy_max, EU(4570.9))
    @test unit_isapprox(
        measurements[end].flux,
        6.15e-10 * inv(one(EU) * NU.m^2 * NU.s),
    )
end
