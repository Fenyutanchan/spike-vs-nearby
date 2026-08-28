# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using QuadGK: quadgk

@testset "continuous PWN physics" begin
    spectrum = InjectionSpectrum_PowerLawWithExponentialCutoff(
        1.5,
        EU(1.0),
        EU(1.0e4);
        cutoff_index=1.0,
        reference_energy=EU(1.0),
    )
    birth_luminosity = EU(10.0, 2)
    spin_down_timescale = EU(100.0, -1)
    efficiency = 0.1
    distance = EU(2.0, -1)
    source_age = EU(10.0, -1)
    source = PWNSource(
        spectrum,
        birth_luminosity,
        spin_down_timescale,
        efficiency,
        distance,
        source_age,
    )

    @testset "pure-dipole spin-down luminosity" begin
        @test unit_isapprox(
            pwn_spin_down_luminosity(source, EU(0.0, -1)),
            birth_luminosity,
        )
        @test unit_isapprox(
            pwn_spin_down_luminosity(source, spin_down_timescale),
            birth_luminosity / 4,
        )

        present_luminosity = EU(2.0, 2)
        pulsar_age = EU(30.0, -1)
        from_present = PWNSource(
            spectrum,
            efficiency,
            distance,
            pulsar_age;
            present_luminosity=present_luminosity,
            spin_down_timescale=EU(10.0, -1),
        )
        @test unit_isapprox(
            pwn_spin_down_luminosity(from_present, pulsar_age),
            present_luminosity,
        )
    end

    @testset "injection power and support" begin
        injection_time = EU(5.0, -1)
        minimum_energy = minimum_injection_energy(spectrum)
        injection_power_value, _ = quadgk(
            EUval(EU, minimum_energy),
            Inf;
            rtol=1e-10,
        ) do energy
            energy * pwn_injection_rate(source, EU(energy), injection_time)
        end
        injection_power = EU(injection_power_value, 2)
        @test unit_isapprox(
            injection_power,
            efficiency * pwn_spin_down_luminosity(source, injection_time);
            rtol=1e-9,
        )
        @test pwn_injection_rate(
            source,
            minimum_energy / 2,
            injection_time,
        ) == 0
        @test pwn_injection_rate(
            source,
            EU(10.0),
            source_age + EU(1.0, -1),
        ) == 0
    end

    @testset "propagated number density" begin
        transport = TransportModel_Homogeneous(
            DiffusionModel_PowerLaw(EU(1.0, -1), 0.0, EU(1.0)),
            EnergyLossModel_PowerLaw(EU(1.0e-3, 2), 2.0, EU(1.0)),
        )
        observed = EU(10.0)
        density = pwn_number_density(source, transport, observed; rtol=1e-10)
        @test EUval(EU, density) > 0

        source_maximum = source_energy_ceiling(transport, observed, source_age)
        direct_integral, _ = quadgk(
            EUval(EU, observed),
            EUval(EU, source_maximum);
            rtol=1e-10,
        ) do source_energy_value
            source_energy = EU(source_energy_value)
            retarded_time = source_age -
                            cooling_time(transport, source_energy, observed)
            pwn_injection_rate(source, source_energy, retarded_time) *
            EUval(
                EU,
                green_kernel(
                    diffusion_loss_integral(
                        transport,
                        source_energy,
                        observed,
                    ),
                    distance,
                ),
            )
        end
        expected_density = EU(direct_integral, 4) /
                           energy_loss_rate(transport, observed)
        @test unit_isapprox(density, expected_density; rtol=1e-10)
    end
end
