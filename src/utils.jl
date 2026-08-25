# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

##############################################################################
function _require_finite(value::Real, name::AbstractString)
    isfinite(value) || throw(DomainError(value, "$name must be finite"))
    return value
end

function _require_finite(value::EnergyUnit, name::AbstractString)
    isfinite(EUval(value)) || throw(DomainError(value, "$name must be finite"))
    return value
end

function _require_positive_finite(value::Real, name::AbstractString)
    _require_finite(value, name)
    value > 0 || throw(DomainError(value, "$name must be positive"))
    return value
end

function _require_positive_finite(value::EnergyUnit, name::AbstractString)
    _require_finite(value, name)
    EUval(value) > 0 || throw(DomainError(value, "$name must be positive"))
    return value
end

function _require_nonnegative_finite(value::Real, name::AbstractString)
    _require_finite(value, name)
    value >= 0 || throw(DomainError(value, "$name must be nonnegative"))
    return value
end

function _require_dimension(
    value::EnergyUnit,
    expected_dimension::Real,
    name::AbstractString,
)
    EUdim(value) == expected_dimension || throw(DimensionMismatch(
        "$name must have natural-unit mass dimension $expected_dimension; " *
        "got $(EUdim(value))",
    ))
    return value
end

function _canonical_unit(
    value::EnergyUnit,
    dimension::Real,
    name::AbstractString,
)
    _require_dimension(value, dimension, name)
    return convert(EU, value)
end
##############################################################################
