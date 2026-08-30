# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

##############################################################################
const _DAMPE_ELECTRON_POSITRON_FLUX_FILE =
    joinpath(external_data_directory, "dampe_e±.csv")

const _DAMPE_ELECTRON_POSITRON_HEADER = [
    "energy_min GeV",
    "energy_max GeV",
    "energy_mean GeV",
    "energy_mean_error GeV",
    "acceptance m^2sr",
    "acceptance_error m^2sr",
    "number_of_events",
    "background_fraction percent",
    "background_fraction_error percent",
    "electron_positron_flux m^-2sr^-1s^-1GeV^-1",
    "electron_positron_flux_error_statistical m^-2sr^-1s^-1GeV^-1",
    "electron_positron_flux_error_systematic m^-2sr^-1s^-1GeV^-1",
]
##############################################################################


##############################################################################
export read_DAMPE_electron_positron_flux

"""
    read_DAMPE_electron_positron_flux([path])

Read the DAMPE combined ``e^-+e^+`` differential-flux table.

The default input is the local transcription `dampe_e±.csv` of Table 1 in the
2017 DAMPE CRE publication. The returned vector contains all 38 published
bins, including separate statistical and systematic flux uncertainties, with
natural units attached.
"""
function read_DAMPE_electron_positron_flux(
    path::AbstractString=_DAMPE_ELECTRON_POSITRON_FLUX_FILE,
)
    data = _read_numeric_csv(path, _DAMPE_ELECTRON_POSITRON_HEADER)
    return [
        _differential_flux_measurement(row[1:4]..., row[10:12]...)
            for row in eachrow(data)
    ]
end
##############################################################################
