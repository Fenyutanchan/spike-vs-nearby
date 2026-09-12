# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ForwardDiff
using NaturalUnits: GeV, MeV

@testset "Empirical background flux" begin
    electron = (C=200.0 * differential_flux_unit, Ebr=GeV(500.0), γ=3.0, Δγ=1.5)
    positron = (C=4.0 * differential_flux_unit, Ebr=GeV(300.0), γ=2.7, Δγ=2.0)

    @testset "Charge spectra below, at and above the break" begin
        for param ∈ (electron, positron),
                energy ∈ (param.Ebr / 10, param.Ebr, 10param.Ebr)
            expected = param.C * (energy / GeV(1.0))^(-param.γ) /
                       (1 + (energy / param.Ebr)^param.Δγ)
            @test unit_isapprox(background_flux(energy, param), expected)
        end
    end

    @testset "Units and reference energy" begin
        energy, reference = GeV(100.0), GeV(10.0)
        converted = map(value -> value isa EnergyUnit ?
                            convert(MeV, value) : value, electron)
        @test unit_isapprox(background_flux(MeV(100000.0), converted;
                                            reference=MeV(10000.0)),
                            background_flux(energy, electron; reference))
        for param ∈ (electron, positron)
            expected = param.C * (energy / reference)^(-param.γ) /
                       (1 + (energy / param.Ebr)^param.Δγ)
            @test unit_isapprox(background_flux(energy, param; reference),
                                expected)
            rescaled = merge(param,
                                (C=param.C * (GeV(1.0) / reference)^param.γ,))
            @test unit_isapprox(background_flux(energy, rescaled; reference),
                                background_flux(energy, param))
        end
    end

    @testset "Spectral slopes and parameter derivatives for \
                uncertainty propagation" begin
        values = [200.0, 500.0, 3.0, 1.5]
        reference = GeV(10.0)
        for energy_value ∈ (100.0, 500.0, 2000.0)
            ratio = energy_value / values[2]
            weight = ratio^values[4] / (1 + ratio^values[4])
            expected_flux = values[1] * (energy_value / 10.0)^(-values[3]) /
                            (1 + ratio^values[4])
            energy_derivative = ForwardDiff.derivative(energy_value) do value
                flux_value(background_flux(GeV(value), electron; reference))
            end
            @test energy_derivative ≈
                    -expected_flux * (values[3] + values[4] * weight) /
                                                                    energy_value
            gradient = ForwardDiff.gradient(values) do param
                physical = (C=param[1] * differential_flux_unit,
                            Ebr=GeV(param[2]), γ=param[3], Δγ=param[4])
                flux_value(background_flux(GeV(energy_value), physical;
                                                                    reference))
            end
            expected_gradient = expected_flux .* [
                inv(values[1]), values[4] * weight / values[2],
                -log(energy_value / 10.0), -weight * log(ratio),
            ]
            @test gradient ≈ expected_gradient
        end
    end
end
