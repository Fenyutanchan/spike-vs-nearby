# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

##############################################################################
export AbstractDiffusionModel, AbstractEnergyLossModel
export AbstractTransportModel

"""
    AbstractDiffusionModel

Supertype for spatially homogeneous diffusion laws.

Concrete subtypes must implement `diffusion_coefficient(model, energy)` for
an `EnergyUnit` energy with natural-unit mass dimension ``+1`` and return an
`EnergyUnit` diffusion coefficient with mass dimension ``-1``.
"""
abstract type AbstractDiffusionModel end

"""
    AbstractEnergyLossModel

Supertype for continuous energy-loss laws.

Concrete subtypes must implement `energy_loss_rate(model, energy)`, returning
the positive loss magnitude

```math
b(E) = -\\frac{\\mathrm{d}E}{\\mathrm{d}t} > 0.
```

The input energy and result must be `EnergyUnit` values with natural-unit mass
dimensions ``+1`` and ``+2``, respectively.
"""
abstract type AbstractEnergyLossModel end

"""
    AbstractTransportModel

Supertype for complete transport models.

Concrete subtypes determine which diffusion, energy-loss, and propagation
methods are available. In particular, a spatially dependent transport model
need not support the homogeneous Syrovatskii-variable interface.
"""
abstract type AbstractTransportModel end
##############################################################################

##############################################################################
export diffusion_coefficient, energy_loss_rate

"""
    diffusion_coefficient(model, energy)

Return the diffusion coefficient at `energy`.

`energy` must have natural-unit mass dimension ``+1``. The result must be an
`EnergyUnit` with mass dimension ``-1``, corresponding to length squared per
time.
"""
function diffusion_coefficient end

"""
    energy_loss_rate(model, energy)

Return the positive continuous energy-loss magnitude

```math
b(E) = -\\frac{\\mathrm{d}E}{\\mathrm{d}t} > 0.
```

`energy` must have natural-unit mass dimension ``+1``. The result must be an
`EnergyUnit` with mass dimension ``+2``, corresponding to energy per time.
"""
function energy_loss_rate end
##############################################################################

##############################################################################
export cooling_time, maximum_cooling_time, source_energy_ceiling

"""
    cooling_time(model, E_source, E_observed)

Return the cooling time

```math
\\tau(E_\\mathrm{s}, E) = \\int_E^{E_\\mathrm{s}} \\frac{\\mathrm{d}E'}{b(E')}
```

from `E_source` to `E_observed`. Both energies must have natural-unit mass
dimension ``+1``. The result must have mass dimension ``-1``.
"""
function cooling_time end

"""
    maximum_cooling_time(model, E_observed)

Return the cooling time from infinite source energy,

```math
\\tau_\\mathrm{max}(E) = \\int_E^\\infty \\frac{\\mathrm{d}E'}{b(E')}.
```

The result must have natural-unit mass dimension ``-1`` and may be infinite.
"""
function maximum_cooling_time end

"""
    source_energy_ceiling(model, E_observed, elapsed_time)

Return the largest source energy compatible with `elapsed_time`, defined by

```math
\\tau(E_\\mathrm{s}, E) = \\Delta t.
```

`E_observed` and the result have natural-unit mass dimension ``+1``;
`elapsed_time` has mass dimension ``-1``. The result may be infinite.
"""
function source_energy_ceiling end
##############################################################################

##############################################################################
export diffusion_loss_integral, maximum_diffusion_loss_integral

"""
    diffusion_loss_integral(model, E_source, E_observed)

Return the diffusion-loss integral

```math
\\Delta\\widetilde{\\lambda}(E_\\mathrm{s}, E)
= \\int_E^{E_\\mathrm{s}} \\frac{D(E')}{b(E')}\\,\\mathrm{d}E'.
```

Both energies must have natural-unit mass dimension ``+1``. The result must
have mass dimension ``-2``, corresponding to an area.
"""
function diffusion_loss_integral end

"""
    maximum_diffusion_loss_integral(
        model,
        E_observed;
        source_energy_max=EU(Inf),
    )

Return the largest diffusion-loss integral allowed by the source-energy
ceiling. The result must have natural-unit mass dimension ``-2`` and may be
infinite.
"""
function maximum_diffusion_loss_integral end
##############################################################################


##############################################################################
export DiffusionModel_PowerLaw, EnergyLossModel_PowerLaw
export TransportModel_Homogeneous

"""
    DiffusionModel_PowerLaw(normalization, index, reference_energy)

Power-law diffusion model

```math
D(E) = D_0 \\left(\\frac{E}{E_0}\\right)^\\delta.
```

`normalization` is an `EnergyUnit` with natural-unit mass dimension ``-1``,
`index` is ``\\delta``, and `reference_energy` is an `EnergyUnit` with mass
dimension ``+1``. For example, a conventional diffusion normalization can be
written as `3e28 * NU.cm^2 / NU.s`.

All dimensionful fields are stored in the canonical [`EU`](@ref) basis. The
underlying numeric values and `index` are promoted together before conversion
to floating point, so a `BigFloat` input preserves arbitrary precision.
"""
struct DiffusionModel_PowerLaw{T<:Real} <: AbstractDiffusionModel
    normalization::EU
    index::T
    reference_energy::EU

    function DiffusionModel_PowerLaw{T}(
        normalization::EU,
        index::T,
        reference_energy::EU,
    ) where {T<:Real}
        _require_dimension(normalization, -1, "diffusion normalization")
        _require_positive_finite(normalization, "diffusion normalization")
        _require_finite(index, "diffusion index")
        _require_dimension(reference_energy, 1, "diffusion reference energy")
        _require_positive_finite(reference_energy, "diffusion reference energy")
        return new{T}(normalization, index, reference_energy)
    end
end

function DiffusionModel_PowerLaw(
    normalization::EnergyUnit,
    index::Real,
    reference_energy::EnergyUnit,
)
    D0_unit = _canonical_unit(normalization, -1, "diffusion normalization")
    E0_unit = _canonical_unit(reference_energy, 1, "diffusion reference energy")
    _require_positive_finite(D0_unit, "diffusion normalization")
    _require_positive_finite(E0_unit, "diffusion reference energy")
    _require_finite(index, "diffusion index")

    D0, delta, E0 = promote(
        EUval(D0_unit),
        index,
        EUval(E0_unit),
    )
    D0, delta, E0 = float(D0), float(delta), float(E0)
    return DiffusionModel_PowerLaw{typeof(delta)}(
        EU(D0, -1),
        delta,
        EU(E0),
    )
end

"""
    EnergyLossModel_PowerLaw(normalization, index, reference_energy)

Power-law energy-loss model

```math
b(E) = b_0 \\left(\\frac{E}{E_0}\\right)^\\alpha,
\\qquad
b(E) = -\\frac{\\mathrm{d}E}{\\mathrm{d}t} > 0.
```

`normalization` is an `EnergyUnit` with natural-unit mass dimension ``+2``,
`index` is ``\\alpha``, and `reference_energy` is an `EnergyUnit` with mass
dimension ``+1``. For example, a conventional loss normalization can be
written as `1e-16 * EU() / NU.s`. The Thomson-limit inverse-Compton plus
synchrotron baseline uses `index = 2`.

All dimensionful fields are stored in the canonical [`EU`](@ref) basis. The
underlying numeric values and `index` are promoted together before conversion
to floating point, so a `BigFloat` input preserves arbitrary precision.
"""
struct EnergyLossModel_PowerLaw{T<:Real} <: AbstractEnergyLossModel
    normalization::EU
    index::T
    reference_energy::EU

    function EnergyLossModel_PowerLaw{T}(
        normalization::EU,
        index::T,
        reference_energy::EU,
    ) where {T<:Real}
        _require_dimension(normalization, 2, "energy-loss normalization")
        _require_positive_finite(normalization, "energy-loss normalization")
        _require_finite(index, "energy-loss index")
        _require_dimension(reference_energy, 1, "energy-loss reference energy")
        _require_positive_finite(reference_energy, "energy-loss reference energy")
        return new{T}(normalization, index, reference_energy)
    end
end

function EnergyLossModel_PowerLaw(
    normalization::EnergyUnit,
    index::Real,
    reference_energy::EnergyUnit,
)
    b0_unit = _canonical_unit(normalization, 2, "energy-loss normalization")
    E0_unit = _canonical_unit(reference_energy, 1, "energy-loss reference energy")
    _require_positive_finite(b0_unit, "energy-loss normalization")
    _require_positive_finite(E0_unit, "energy-loss reference energy")
    _require_finite(index, "energy-loss index")

    b0, alpha, E0 = promote(
        EUval(b0_unit),
        index,
        EUval(E0_unit),
    )
    b0, alpha, E0 = float(b0), float(alpha), float(E0)
    return EnergyLossModel_PowerLaw{typeof(alpha)}(
        EU(b0, 2),
        alpha,
        EU(E0),
    )
end

"""
    TransportModel_Homogeneous(diffusion, energy_loss)

Combine a homogeneous diffusion law with a continuous energy-loss law.

The Green kernel for this model assumes an infinite homogeneous medium. A
spatially dependent or two-zone diffusion coefficient requires a different
transport solver.
"""
struct TransportModel_Homogeneous{
    D<:AbstractDiffusionModel,
    L<:AbstractEnergyLossModel,
} <: AbstractTransportModel
    diffusion::D
    energy_loss::L
end
##############################################################################


##############################################################################
"""
    diffusion_coefficient(
        model::DiffusionModel_PowerLaw,
        energy::EnergyUnit,
    )

Evaluate the power-law diffusion coefficient at `energy`.

`energy` must have natural-unit mass dimension ``+1``. The result is an
`EnergyUnit` with mass dimension ``-1`` (equivalent to length squared per
time) in the canonical [`EU`](@ref) basis.
"""
function diffusion_coefficient(
    model::DiffusionModel_PowerLaw,
    energy::EnergyUnit,
)
    energy = _canonical_unit(energy, 1, "energy")
    _require_positive_finite(energy, "energy")
    return model.normalization *
           (energy / model.reference_energy)^model.index
end

diffusion_coefficient(model::TransportModel_Homogeneous, energy::EnergyUnit) =
    diffusion_coefficient(model.diffusion, energy)

"""
    energy_loss_rate(
        model::EnergyLossModel_PowerLaw,
        energy::EnergyUnit,
    )

Evaluate the positive power-law loss magnitude

```math
b(E) = -\\frac{\\mathrm{d}E}{\\mathrm{d}t}
```

at `energy`. The input energy must have natural-unit mass dimension ``+1``;
the result has mass dimension ``+2`` (equivalent to energy per time) in the
canonical [`EU`](@ref) basis.
"""
function energy_loss_rate(
    model::EnergyLossModel_PowerLaw,
    energy::EnergyUnit,
)
    energy = _canonical_unit(energy, 1, "energy")
    _require_positive_finite(energy, "energy")
    return model.normalization *
           (energy / model.reference_energy)^model.index
end

energy_loss_rate(model::TransportModel_Homogeneous, energy::EnergyUnit) =
    energy_loss_rate(model.energy_loss, energy)
##############################################################################


##############################################################################
"""
    cooling_time(
        model::TransportModel_Homogeneous,
        E_source::EnergyUnit,
        E_observed::EnergyUnit,
    )

Evaluate the defining cooling-time integral numerically with `quadgk`:

```math
\\tau(E_\\mathrm{s}, E)
    = \\int_E^{E_\\mathrm{s}} \\frac{\\mathrm{d}E'}{b(E')}
```

This fallback applies to homogeneous energy-loss models without a specialized
analytic method.
"""
function cooling_time(
    model::TransportModel_Homogeneous,
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
)
    source, observed = _check_energy_interval(E_source, E_observed)
    source == observed && return _zero_unit_like(observed, -1)

    value = _integrate_over_energy(
        energy -> inv(energy_loss_rate(model, energy)),
        source,
        observed,
        -1,
    )
    return _require_nonnegative_finite(value, "cooling time")
end

"""
    cooling_time(
        model::TransportModel_Homogeneous{D,L},
        E_source::EnergyUnit,
        E_observed::EnergyUnit,
    ) where {D<:AbstractDiffusionModel,L<:EnergyLossModel_PowerLaw}

Evaluate the cooling time analytically for a power-law energy-loss model. The
implementation uses `log1p` and `_exprel` to remain stable when the two
energies are close and when the power-law integral becomes logarithmic.
"""
function cooling_time(
    model::TransportModel_Homogeneous{D, L},
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
) where {D<:AbstractDiffusionModel, L<:EnergyLossModel_PowerLaw}
    source, observed = _check_energy_interval(E_source, E_observed)
    log_ratio = _energy_log_ratio(source, observed)
    alpha = model.energy_loss.index
    value = observed / energy_loss_rate(model, observed) *
            log_ratio * _exprel((1 - alpha) * log_ratio)
    return _require_nonnegative_finite(value, "cooling time")
end

"""
    maximum_cooling_time(
        model::TransportModel_Homogeneous{D,L},
        E_observed::EnergyUnit,
    ) where {D<:AbstractDiffusionModel,L<:EnergyLossModel_PowerLaw}

Evaluate

```math
\\tau_\\mathrm{max}(E)
= \\int_E^\\infty \\frac{\\mathrm{d}E'}{b(E')}
```

The returned `EnergyUnit` has natural-unit mass dimension ``-1``. The result
is finite only for a power-law loss index greater than one.
"""
function maximum_cooling_time(
    model::TransportModel_Homogeneous{D, L},
    E_observed::EnergyUnit,
) where {D<:AbstractDiffusionModel, L<:EnergyLossModel_PowerLaw}
    observed = _canonical_unit(E_observed, 1, "observed energy")
    _require_positive_finite(observed, "observed energy")
    alpha = model.energy_loss.index
    alpha > 1 || return _infinite_unit_like(observed, -1)
    return observed / energy_loss_rate(model, observed) / (alpha - 1)
end

"""
    source_energy_ceiling(
        model::TransportModel_Homogeneous{D,L},
        E_observed::EnergyUnit,
        elapsed_time::EnergyUnit,
    ) where {D<:AbstractDiffusionModel,L<:EnergyLossModel_PowerLaw}

Invert the power-law cooling-time relation analytically. The result solves

```math
\\tau(E_\\mathrm{s}, E) = \\Delta t.
```

An infinite source energy is returned when `elapsed_time` reaches or exceeds
the cooling time from infinite source energy.
"""
function source_energy_ceiling(
    model::TransportModel_Homogeneous{D, L},
    E_observed::EnergyUnit,
    elapsed_time::EnergyUnit,
) where {D<:AbstractDiffusionModel, L<:EnergyLossModel_PowerLaw}
    observed = _canonical_unit(E_observed, 1, "observed energy")
    elapsed = _canonical_unit(elapsed_time, -1, "elapsed time")
    _require_positive_finite(observed, "observed energy")
    _require_nonnegative_finite(elapsed, "elapsed time")
    iszero(elapsed) && return observed

    maximum_time = maximum_cooling_time(model, observed)
    elapsed >= maximum_time && return _infinite_unit_like(observed, 1)

    alpha = model.energy_loss.index
    scaled_time = energy_loss_rate(model, observed) / observed * elapsed
    log_ratio = if alpha == 1
        scaled_time
    else
        log1p((1 - alpha) * scaled_time) / (1 - alpha)
    end

    source_value = EUval(observed) * exp(log_ratio)
    isfinite(source_value) || return _infinite_unit_like(observed, 1)
    return EU(source_value)
end
##############################################################################


##############################################################################
"""
    diffusion_loss_integral(
        model::TransportModel_Homogeneous,
        E_source::EnergyUnit,
        E_observed::EnergyUnit,
    )

Evaluate the defining diffusion-loss integral numerically with `quadgk`:

```math
\\Delta\\widetilde{\\lambda}(E_\\mathrm{s}, E)
= \\int_E^{E_\\mathrm{s}} \\frac{D(E')}{b(E')}\\,\\mathrm{d}E'
```

This fallback applies when the diffusion and energy-loss pair has no
specialized analytic method.
"""
function diffusion_loss_integral(
    model::TransportModel_Homogeneous,
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
)
    source, observed = _check_energy_interval(E_source, E_observed)
    source == observed && return _zero_unit_like(observed, -2)

    value = _integrate_over_energy(
        energy -> diffusion_coefficient(model, energy) /
                  energy_loss_rate(model, energy),
        source,
        observed,
        -2,
    )
    return _require_nonnegative_finite(value, "diffusion-loss integral")
end

"""
    diffusion_loss_integral(
        model::TransportModel_Homogeneous{D,L},
        E_source::EnergyUnit,
        E_observed::EnergyUnit,
    ) where {D<:DiffusionModel_PowerLaw,L<:EnergyLossModel_PowerLaw}

Evaluate the diffusion-loss integral analytically for a power-law diffusion
and loss pair. The implementation uses `log1p` and `_exprel` to remain stable
near equal energies and in the logarithmic limit
``1 + \\delta - \\alpha = 0``.
"""
function diffusion_loss_integral(
    model::TransportModel_Homogeneous{D,L},
    E_source::EnergyUnit,
    E_observed::EnergyUnit,
) where {D<:DiffusionModel_PowerLaw,L<:EnergyLossModel_PowerLaw}
    source, observed = _check_energy_interval(E_source, E_observed)
    log_ratio = _energy_log_ratio(source, observed)
    delta = model.diffusion.index
    alpha = model.energy_loss.index
    value = observed * diffusion_coefficient(model, observed) /
            energy_loss_rate(model, observed) * log_ratio *
            _exprel((1 + delta - alpha) * log_ratio)
    return _require_nonnegative_finite(value, "diffusion-loss integral")
end

"""
    maximum_diffusion_loss_integral(
        model::TransportModel_Homogeneous{D,L},
        E_observed::EnergyUnit;
        source_energy_max::EnergyUnit=EU(Inf),
    ) where {D<:DiffusionModel_PowerLaw,L<:EnergyLossModel_PowerLaw}

Return the maximum diffusion-loss integral at `E_observed` as an `EnergyUnit`
with natural-unit mass dimension ``-2``.

For an infinite source-energy ceiling, convergence requires
``\\alpha - \\delta > 1``. Set a finite `source_energy_max` to evaluate the
diffusion-loss integral only up to that injection cutoff.
"""
function maximum_diffusion_loss_integral(
    model::TransportModel_Homogeneous{D,L},
    E_observed::EnergyUnit;
    source_energy_max::EnergyUnit=EU(Inf),
) where {D<:DiffusionModel_PowerLaw,L<:EnergyLossModel_PowerLaw}
    observed = _canonical_unit(E_observed, 1, "observed energy")
    source_maximum = _canonical_unit(
        source_energy_max,
        1,
        "source-energy maximum",
    )
    _require_positive_finite(observed, "observed energy")

    if isinf(source_maximum)
        EUval(source_maximum) > 0 || throw(DomainError(
            source_maximum,
            "source-energy maximum must be finite or +Inf",
        ))
    elseif isnan(source_maximum)
        throw(DomainError(
            source_maximum,
            "source-energy maximum must be finite or +Inf",
        ))
    else
        return diffusion_loss_integral(model, source_maximum, observed)
    end

    convergence_power = model.energy_loss.index - model.diffusion.index - 1
    convergence_power > 0 || return _infinite_unit_like(observed, -2)
    return observed * diffusion_coefficient(model, observed) /
           energy_loss_rate(model, observed) / convergence_power
end
##############################################################################


##############################################################################
export green_kernel

function _log_green_kernel(
    diffusion_integral::EnergyUnit,
    distance::EnergyUnit,
)
    diffusion_integral = _canonical_unit(
        diffusion_integral,
        -2,
        "diffusion-loss integral",
    )
    distance = _canonical_unit(distance, -1, "distance")
    _require_nonnegative_finite(
        diffusion_integral,
        "diffusion-loss integral",
    )
    _require_nonnegative_finite(distance, "distance")
    lambda, radius = promote(EUval(diffusion_integral), EUval(distance))
    lambda, radius = float(lambda), float(radius)

    if iszero(lambda)
        iszero(radius) && throw(DomainError(
            (diffusion_integral, distance),
            "the zero-distance, zero-integral Green kernel is a " *
            "delta-distribution singularity",
        ))
        return oftype(lambda, -Inf)
    end

    log_four_pi = log(oftype(lambda, 4 * pi))
    scaled_radius = radius / (2 * sqrt(lambda))
    return oftype(lambda, -3 / 2) * (log_four_pi + log(lambda)) -
           scaled_radius^2
end

"""
    green_kernel(diffusion_integral, distance)

Return the normalized three-dimensional Green kernel as an `EnergyUnit` with
natural-unit mass dimension ``+3`` (an inverse volume).
"""
green_kernel(diffusion_integral::EnergyUnit, distance::EnergyUnit) =
    EU(exp(_log_green_kernel(diffusion_integral, distance)), 3)
##############################################################################
