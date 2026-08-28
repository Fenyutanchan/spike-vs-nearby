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
