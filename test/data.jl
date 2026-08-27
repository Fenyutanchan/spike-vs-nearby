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

    frequency = parse(Float64, Geminga_ATNF["F0"].value)
    frequency_derivative = parse(Float64, Geminga_ATNF["F1"].value)
    characteristic_age = -frequency / (2 * frequency_derivative)
    @test characteristic_age / (365.25 * 86400) ≈ 342_000 rtol = 2e-3

    pulsar = read_ATNF_pulsar("B0011+47")
    @test pulsar.name == "J0014+4746"
    @test pulsar["PSRB"].value == "B0011+47"
end
