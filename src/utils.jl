# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using QuadGK: quadgk

##############################################################################
@inline function _energy_log_ratio(
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
)
    return log1p((E_source - E_observed) / E_observed)
end

@inline function _exprel(value::Real)
    abs(value) < 1e-5 || return expm1(value) / value
    return evalpoly(
        value,
        (one(value), one(value) / 2, one(value) / 6,
         one(value) / 24, one(value) / 120, one(value) / 720),
    )
end

function _check_energy_interval(
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
)
    source = _canonical_unit(E_source, 1, "source energy")
    observed = _canonical_unit(E_observed, 1, "observed energy")
    _require_positive_finite(source, "source energy")
    _require_positive_finite(observed, "observed energy")
    source >= observed || throw(DomainError(
        (source, observed),
        "source energy must be greater than or equal to observed energy",
    ))
    return source, observed
end

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

function _require_nonnegative_finite(value::EnergyUnit, name::AbstractString)
    _require_finite(value, name)
    EUval(value) >= 0 || throw(DomainError(value, "$name must be nonnegative"))
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

@inline function _zero_unit_like(reference::EnergyUnit, dimension::Real)
    value = zero(float(EUval(EU, reference)))
    return EU(value, dimension)
end

@inline function _infinite_unit_like(reference::EnergyUnit, dimension::Real)
    value = oftype(float(EUval(EU, reference)), Inf)
    return EU(value, dimension)
end

function _integrate_over_energy(
    integrand,
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
    result_dimension::Real,
)
    source, observed = _check_energy_interval(E_source, E_observed)
    source_value, observed_value = promote(
        EUval(EU, source),
        EUval(EU, observed),
    )
    source_value = float(source_value)
    observed_value = float(observed_value)
    integrand_dimension = result_dimension - 1

    value, _ = quadgk(observed_value, source_value) do energy_value
        integrand_value = integrand(EU(energy_value))
        _require_dimension(
            integrand_value,
            integrand_dimension,
            "energy integrand",
        )
        return EUval(EU, integrand_value)
    end
    return EU(value, result_dimension)
end
##############################################################################
