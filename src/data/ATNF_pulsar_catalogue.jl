# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Artifacts

################################################################################
export ATNF_PULSAR_CATALOGUE_DATABASE

const _ATNF_PULSAR_CATALOGUE_VERSION = v"2.8.1"
const _ATNF_PULSAR_CATALOGUE_ARTIFACT_NAME =
    "ATNF_Pulsar_Catalogue_v$(_ATNF_PULSAR_CATALOGUE_VERSION)"

"""
    ATNF_PULSAR_CATALOGUE_DATABASE

Path to the `psrcat.db` file in the pinned ATNF catalogue artifact.
"""
const ATNF_PULSAR_CATALOGUE_DATABASE = let
    artifacts_toml = joinpath(project_directory, "Artifacts.toml")
    hash = artifact_hash(
        _ATNF_PULSAR_CATALOGUE_ARTIFACT_NAME,
        artifacts_toml,
    )
    joinpath(artifact_path(hash), "psrcat_tar", "psrcat.db")
end
################################################################################

################################################################################
export ATNFParameter, ATNFPulsar

"""
    ATNFParameter

One parameter read from an ATNF `psrcat.db` record.

`value` preserves the catalogue token exactly, and `uncertainty` preserves the
optional uncertainty in the last quoted digits. `unit` records the display
unit specified by `psrcat` when recognized by this lightweight parser, while
`reference` preserves the optional ATNF bibliography key. No derived
quantities are assigned by this layer.
"""
struct ATNFParameter
    value::String
    uncertainty::Union{Nothing, String}
    unit::Union{Nothing, String}
    reference::Union{Nothing, String}
end

"""
    ATNFPulsar

Lightweight representation of one pulsar record in `psrcat.db`.

`name` is the `PSRJ` value when available and otherwise the `PSRB` value.
Parameters are indexed case-insensitively by their ATNF labels, for example
`pulsar["F0"]` or `pulsar["dist_a"]`.
"""
struct ATNFPulsar
    name::String
    parameters::Dict{String, ATNFParameter}
end

Base.getindex(pulsar::ATNFPulsar, parameter::AbstractString) =
    pulsar.parameters[uppercase(parameter)]

Base.haskey(pulsar::ATNFPulsar, parameter::AbstractString) =
    haskey(pulsar.parameters, uppercase(parameter))

function Base.show(io::IO, pulsar::ATNFPulsar)
    print(io,
        "ATNFPulsar(",
        repr(pulsar.name), ", ",
        length(pulsar.parameters),
        " parameters)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", pulsar::ATNFPulsar)
    labels = sort!(collect(keys(pulsar.parameters)))
    label_width = maximum(textwidth, labels; init=0)

    print(io, "ATNFPulsar ", pulsar.name, " (", length(labels), " parameters)")

    for label in labels
        parameter = pulsar.parameters[label]
        measurement = isnothing(parameter.uncertainty) ?
            parameter.value : "$(parameter.value)($(parameter.uncertainty))"

        print(io, '\n', "  ", rpad(label, label_width + 2), measurement)
        isnothing(parameter.unit) || print(io, ' ', parameter.unit)
        isnothing(parameter.reference) ||
            print(io, " [", parameter.reference, ']')
    end
end
################################################################################


################################################################################
export read_ATNF_pulsar

"""
    read_ATNF_pulsar(name; database=ATNF_PULSAR_CATALOGUE_DATABASE)

Read one pulsar from an ATNF `psrcat.db` file by its `PSRJ` or `PSRB` name.

The lookup is case-insensitive. The parser retains the raw value, uncertainty,
and reference tokens supplied by ATNF, ignores comments, and follows `psrcat`
in keeping the last occurrence of a repeated parameter. It does not reproduce
the derived-parameter machinery of the `psrcat` executable.

Raise `KeyError` when neither catalogue name matches `name`.
"""
function read_ATNF_pulsar(
    name::AbstractString;
    database::AbstractString=ATNF_PULSAR_CATALOGUE_DATABASE,
)
    target = uppercase(strip(name))
    parameters = Dict{String,ATNFParameter}()

    pulsar = open(database, "r") do input
        for line in eachline(input)
            stripped = strip(line)

            if startswith(stripped, '@')
                pulsar = _matching_ATNF_pulsar(parameters, target)
                isnothing(pulsar) || return pulsar
                empty!(parameters)
            elseif isempty(stripped) || startswith(stripped, '#')
                continue
            else
                parameter = _parse_ATNF_parameter(stripped)
                isnothing(parameter) && continue
                label, entry = parameter
                parameters[label] = entry
            end
        end

        return _matching_ATNF_pulsar(parameters, target)
    end

    isnothing(pulsar) || return pulsar
    throw(KeyError(name))
end

function _parse_ATNF_parameter(line::AbstractString)
    columns = split(line)
    length(columns) >= 2 || return nothing

    label = uppercase(columns[1])
    value = columns[2]
    uppercase(value) == "NULL" && return nothing

    uncertainty = nothing
    reference = nothing
    if length(columns) >= 3
        if _is_ATNF_uncertainty(columns[3])
            uncertainty = String(columns[3])
            length(columns) >= 4 && (reference = String(columns[4]))
        else
            reference = String(columns[3])
        end
    end

    return String(label) => ATNFParameter(
        String(value),
        uncertainty,
        _ATNF_parameter_unit(label),
        reference,
    )
end

function _matching_ATNF_pulsar(
    parameters::Dict{String,ATNFParameter},
    target::AbstractString,
)
    isempty(parameters) && return nothing

    for label in ("PSRJ", "PSRB")
        haskey(parameters, label) || continue
        uppercase(parameters[label].value) == target || continue
        canonical_name = haskey(parameters, "PSRJ") ?
                         parameters["PSRJ"].value :
                         parameters["PSRB"].value
        return ATNFPulsar(canonical_name, copy(parameters))
    end

    return nothing
end

_is_ATNF_uncertainty(token::AbstractString) =
    !isempty(token) && isdigit(first(token))

const _ATNF_PARAMETER_UNITS = Dict(
    "A1" => "lt-s",
    "DECJ" => "dms",
    "DM" => "cm^-3 pc",
    "DMEPOCH" => "MJD",
    "ELAT" => "deg",
    "ELONG" => "deg",
    "FINISH" => "MJD",
    "GAMMA" => "s",
    "GB" => "deg",
    "GL" => "deg",
    "OM" => "deg",
    "OMDOT" => "deg/yr",
    "P0" => "s",
    "PEPOCH" => "MJD",
    "PMB" => "mas/yr",
    "PMDEC" => "mas/yr",
    "PMELAT" => "mas/yr",
    "PMELONG" => "mas/yr",
    "PML" => "mas/yr",
    "PMRA" => "mas/yr",
    "POSEPOCH" => "MJD",
    "PX" => "mas",
    "RAJ" => "hms",
    "RM" => "rad m^-2",
    "START" => "MJD",
    "T0" => "MJD",
    "TASC" => "MJD",
    "TAU_SC" => "s",
    "TRES" => "us",
    "TZRFRQ" => "MHz",
    "TZRMJD" => "MJD",
    "W10" => "ms",
    "W50" => "ms",
)

function _ATNF_parameter_unit(label::AbstractString)
    startswith(label, "DIST") && return "kpc"
    occursin(r"^S(?:\d+|\d+G)$", label) && return "mJy"
    occursin(r"^PB(?:_\d+)?$", label) && return "d"
    occursin(r"^A1(?:_\d+)?$", label) && return "lt-s"
    occursin(r"^(?:T0|TASC)(?:_\d+)?$", label) && return "MJD"
    occursin(r"^OM(?:_\d+)?$", label) && return "deg"
    occursin(r"^OMDOT(?:_\d+)?$", label) && return "deg/yr"

    frequency_match = match(r"^F(\d+)$", label)
    if !isnothing(frequency_match)
        order = parse(Int, only(frequency_match.captures))
        return iszero(order) ? "Hz" : "s^-$((order + 1))"
    end

    dispersion_match = match(r"^DM(\d+)$", label)
    if !isnothing(dispersion_match)
        order = parse(Int, only(dispersion_match.captures))
        return "cm^-3 pc yr^-$order"
    end

    return get(_ATNF_PARAMETER_UNITS, label, nothing)
end
################################################################################

################################################################################
export pulsar_characteristic_age, pulsar_spin_down_luminosity

const _FIDUCIAL_NEUTRON_STAR_MOMENT_OF_INERTIA =
    1e45 * NU.g * NU.cm^2

"""
    pulsar_characteristic_age(pulsar)

Return the characteristic age inferred from the central ATNF values of
`F0` and `F1`,

```math
\\tau_\\mathrm{c} = -\\frac{\\nu}{2\\dot{\\nu}}.
```

This estimate assumes magnetic-dipole braking with braking index ``n=3`` and
an initial spin period much shorter than the present period. It is therefore
a model-dependent age proxy rather than a direct age measurement. The result
has natural-unit mass dimension ``-1``.
"""
function pulsar_characteristic_age(pulsar::ATNFPulsar)
    frequency = parse(Float64, pulsar["F0"].value) / NU.s
    frequency_derivative =
        parse(Float64, pulsar["F1"].value) / NU.s^2
    age = -frequency / (2 * frequency_derivative)
    return _require_positive_finite(age, "pulsar characteristic age")
end

"""
    pulsar_spin_down_luminosity(
        pulsar;
        moment_of_inertia=1e45 * NU.g * NU.cm^2,
    )

Return the present spin-down luminosity inferred from the central ATNF values
of `F0` and `F1`,

```math
L_\\mathrm{sd} = -4\\pi^2 I\\nu\\dot{\\nu}.
```

`moment_of_inertia` is the assumed neutron-star moment of inertia ``I`` and
must have natural-unit mass dimension ``-1``. Its default value is the
fiducial ATNF convention ``10^{45}\\,\\mathrm{g\\,cm^2}``. The result has
mass dimension ``+2``.
"""
function pulsar_spin_down_luminosity(
    pulsar::ATNFPulsar;
    moment_of_inertia::EnergyUnit=
        _FIDUCIAL_NEUTRON_STAR_MOMENT_OF_INERTIA,
)
    inertia = _canonical_unit(
        moment_of_inertia,
        -1,
        "neutron-star moment of inertia",
    )
    _require_positive_finite(inertia, "neutron-star moment of inertia")

    frequency = parse(Float64, pulsar["F0"].value) / NU.s
    frequency_derivative =
        parse(Float64, pulsar["F1"].value) / NU.s^2
    luminosity = -4 * pi^2 * inertia * frequency * frequency_derivative
    return _require_positive_finite(luminosity, "pulsar spin-down luminosity")
end
################################################################################

################################################################################
export Geminga_ATNF, Monogem_ATNF

const Geminga_ATNF = read_ATNF_pulsar("J0633+1746")
const Monogem_ATNF = read_ATNF_pulsar("J0659+1414")
################################################################################
