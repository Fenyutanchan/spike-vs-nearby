# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using DelimitedFiles: readdlm

################################################################################
export DifferentialFluxMeasurement

"""
    DifferentialFluxMeasurement

One binned measurement of a differential particle flux.

`energy_min`, `energy_max`, `energy`, and `energy_error` have natural-unit
mass dimension ``+1``. `flux`, `statistical_error`, and `systematic_error`
have mass dimension ``+2``, corresponding to
``\\mathrm{GeV^{-1}\\,m^{-2}\\,s^{-1}\\,sr^{-1}}`` in the experimental tables.
The `energy` field stores the collaboration's representative energy for the
bin rather than an arithmetic bin centre.
"""
struct DifferentialFluxMeasurement
    energy_min::EU
    energy_max::EU
    energy::EU
    energy_error::EU
    flux::EU
    statistical_error::EU
    systematic_error::EU

    function DifferentialFluxMeasurement(
        energy_min::EnergyUnit,
        energy_max::EnergyUnit,
        energy::EnergyUnit,
        energy_error::EnergyUnit,
        flux::EnergyUnit,
        statistical_error::EnergyUnit,
        systematic_error::EnergyUnit,
    )
        lower = _canonical_unit(energy_min, 1, "lower bin energy")
        upper = _canonical_unit(energy_max, 1, "upper bin energy")
        representative = _canonical_unit(energy, 1, "representative energy")
        energy_uncertainty = _canonical_unit(
            energy_error,
            1,
            "representative-energy uncertainty",
        )
        measured_flux = _canonical_unit(flux, 2, "differential flux")
        statistical = _canonical_unit(
            statistical_error,
            2,
            "statistical flux uncertainty",
        )
        systematic = _canonical_unit(
            systematic_error,
            2,
            "systematic flux uncertainty",
        )

        _require_positive_finite(lower, "lower bin energy")
        _require_positive_finite(upper, "upper bin energy")
        lower < representative < upper || throw(DomainError(
            representative,
            "representative energy must lie inside its bin",
        ))
        _require_nonnegative_finite(
            energy_uncertainty,
            "representative-energy uncertainty",
        )
        _require_nonnegative_finite(measured_flux, "differential flux")
        _require_nonnegative_finite(
            statistical,
            "statistical flux uncertainty",
        )
        _require_nonnegative_finite(
            systematic,
            "systematic flux uncertainty",
        )

        return new(
            lower,
            upper,
            representative,
            energy_uncertainty,
            measured_flux,
            statistical,
            systematic,
        )
    end
end

const _EXPERIMENTAL_DIFFERENTIAL_FLUX_UNIT =
    inv(one(EU) * NU.m^2 * NU.s)

function _read_numeric_csv(
    path::AbstractString,
    expected_header::AbstractVector{<:AbstractString},
)
    raw, header = readdlm(path, ','; header=true)
    String.(vec(header)) == expected_header || throw(ArgumentError(
        "unexpected CSV header in $(abspath(path))",
    ))
    return Float64.(raw)
end

function _differential_flux_measurement(
    energy_min::Real,
    energy_max::Real,
    energy::Real,
    energy_error::Real,
    flux::Real,
    statistical_error::Real,
    systematic_error::Real,
)
    return DifferentialFluxMeasurement(
        EU(energy_min),
        EU(energy_max),
        EU(energy),
        EU(energy_error),
        flux * _EXPERIMENTAL_DIFFERENTIAL_FLUX_UNIT,
        statistical_error * _EXPERIMENTAL_DIFFERENTIAL_FLUX_UNIT,
        systematic_error * _EXPERIMENTAL_DIFFERENTIAL_FLUX_UNIT,
    )
end
################################################################################
