# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

################################################################################
const _AMS02_ELECTRON_FLUX_FILE =
    joinpath(external_data_directory, "ams02_e-.csv")
const _AMS02_POSITRON_FLUX_FILE =
    joinpath(external_data_directory, "ams02_e+.csv")

const _AMS02_ELECTRON_HEADER = [
    "energy_min GeV",
    "energy_max GeV",
    "energy_mean GeV",
    "energy_mean_error_systematic_total GeV",
    "number_of_events",
    "number_of_events_error_statistical",
    "electron_flux m^-2sr^-1s^-1GeV^-1",
    "electron_flux_error_statistical m^-2sr^-1s^-1GeV^-1",
    "electron_flux_error_systematic_total m^-2sr^-1s^-1GeV^-1",
]

const _AMS02_POSITRON_HEADER = [
    "energy_min GeV",
    "energy_max GeV",
    "energy_mean GeV",
    "energy_mean_error_systematic_total GeV",
    "number_of_events",
    "number_of_events_error_statistical",
    "positron_flux m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_statistical m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_template_definition m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_charge_confusion m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_efficiency_correction m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_unfolding m^-2sr^-1s^-1GeV^-1",
    "positron_flux_error_systematic_total m^-2sr^-1s^-1GeV^-1",
]
################################################################################


################################################################################
export read_AMS02_electron_flux, read_AMS02_positron_flux

"""
    read_AMS02_electron_flux()

Read the AMS-02 electron differential-flux table.

The input is the collaboration-provided `ams02_e-.csv` stored in
`data/ext`. The returned vector contains the published representative energy,
flux, and separate statistical and total systematic uncertainties with
natural units attached.
"""
function read_AMS02_electron_flux()
    data = _read_numeric_csv(_AMS02_ELECTRON_FLUX_FILE, _AMS02_ELECTRON_HEADER)
    return [
        _differential_flux_measurement(row[1:4]..., row[7:9]...)
            for row in eachrow(data)
    ]
end

"""
    read_AMS02_positron_flux()

Read the AMS-02 positron differential-flux table.

The input is the collaboration-provided `ams02_e+.csv` stored in
`data/ext`. The returned vector retains the published statistical uncertainty
and the collaboration's total systematic uncertainty; its individual
systematic components are not recombined by this lightweight interface.
"""
function read_AMS02_positron_flux()
    data = _read_numeric_csv(_AMS02_POSITRON_FLUX_FILE, _AMS02_POSITRON_HEADER)
    return [
        _differential_flux_measurement(row[1:4]..., row[7], row[8], row[13])
            for row in eachrow(data)
    ]
end
################################################################################


################################################################################
export combine_AMS02_electron_positron_flux

"""
    combine_AMS02_electron_positron_flux()

Read the fixed AMS-02 electron and positron tables and combine matching bins
into the measured ``e^-+e^+`` differential flux.

Statistical and total systematic uncertainties are each combined in quadrature.
The electron and positron errors are assumed to be uncorrelated in each bin.

The AMS-02 inputs share 74 bins through ``1\\,\\mathrm{TeV}``; the final
electron-only bin is deliberately excluded from the combined result.
"""
function combine_AMS02_electron_positron_flux()
    electrons = read_AMS02_electron_flux()
    positrons = read_AMS02_positron_flux()

    number_of_common_bins = min(length(electrons), length(positrons))
    combined = Vector{DifferentialFluxMeasurement}(
        undef,
        number_of_common_bins,
    )

    for index ∈ eachindex(combined)
        electron = electrons[index]
        positron = positrons[index]
        _matching_AMS02_energy_bin(electron, positron)

        combined[index] = DifferentialFluxMeasurement(
            electron.energy_min,
            electron.energy_max,
            electron.energy,
            electron.energy_error,
            electron.flux + positron.flux,
            sqrt(electron.statistical_error^2 + positron.statistical_error^2),
            sqrt(electron.systematic_error^2 + positron.systematic_error^2),
        )
    end

    return combined
end

function _matching_AMS02_energy_bin(
    electron::DifferentialFluxMeasurement,
    positron::DifferentialFluxMeasurement,
)
    electron.energy_min == positron.energy_min &&
    electron.energy_max == positron.energy_max &&
    electron.energy == positron.energy &&
    electron.energy_error == positron.energy_error && return nothing

    ArgumentError(
        "AMS-02 electron and positron energy bins do not match"
    ) |> throw
end

################################################################################
