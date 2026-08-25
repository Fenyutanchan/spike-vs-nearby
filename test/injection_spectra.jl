# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using QuadGK: quadgk

function numerical_injection_energy_integral(
    spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
    ;
    kwargs...,
)
    gamma = spectrum.index
    beta = spectrum.cutoff_index
    lower = spectrum.minimum_energy / spectrum.cutoff_energy
    upper = oftype(float(lower), Inf)

    dimensionless_integral, _ = quadgk(lower, upper; kwargs...) do ratio
        cutoff_power = ratio^beta
        isinf(cutoff_power) && return zero(ratio + gamma + beta)
        return exp((1 - gamma) * log(ratio) - cutoff_power)
    end

    prefactor = spectrum.cutoff_energy^2 *
                (spectrum.cutoff_energy /
                 spectrum.reference_energy)^(-gamma)
    return prefactor * dimensionless_integral
end

@testset "power-law injection spectrum with exponential cutoff" begin
    spectrum = InjectionSpectrum_PowerLawWithExponentialCutoff(
        1.5,
        EU(2.0),
        EU(10.0);
        cutoff_index=2.0,
        reference_energy=EU(3.0),
    )

    @test unit_isapprox(minimum_injection_energy(spectrum), EU(2.0))
    @test injection_shape(spectrum, EU(1.0)) == 0.0
    @test injection_shape(spectrum, EU(2.0)) ≈
          (2 / 3)^(-1.5) * exp(-(2 / 10)^2)

    elementary_spectrum = InjectionSpectrum_PowerLawWithExponentialCutoff(
        0.0,
        EU(2.0),
        EU(10.0),
    )
    normalization = injection_energy_integral(elementary_spectrum)
    expected = EU(10.0^2 * (1 + 2 / 10) * exp(-2 / 10), 2)
    @test unit_isapprox(normalization, expected; rtol=1e-12)

    @testset "incomplete-Gamma result against QuadGK" begin
        cases = (
            (1.8, 1.0, 100.0, 1.0, 1.0),
            (2.0, 3.0, 80.0, 2.0, 5.0),
            (2.4, 10.0, 300.0, 0.75, 2.0),
        )

        for (index, minimum, cutoff, cutoff_index, reference) in cases
            spectrum_case = InjectionSpectrum_PowerLawWithExponentialCutoff(
                index,
                EU(minimum),
                EU(cutoff);
                cutoff_index=cutoff_index,
                reference_energy=EU(reference),
            )
            analytic = injection_energy_integral(spectrum_case)
            numerical = numerical_injection_energy_integral(
                spectrum_case;
                rtol=1e-12,
            )

            @test unit_isapprox(numerical, analytic; rtol=1e-11)
        end
    end
end
