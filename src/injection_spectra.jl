# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using SpecialFunctions: gamma

##############################################################################
export AbstractInjectionSpectrum

"""
    AbstractInjectionSpectrum

Supertype for dimensionless injection-spectrum shapes.

Concrete subtypes must implement `injection_shape(spectrum, energy)`,
`minimum_injection_energy(spectrum)`, and
`injection_energy_integral(spectrum)`. Energies have natural-unit mass
dimension ``+1``, while the normalization integral has mass dimension ``+2``.
"""
abstract type AbstractInjectionSpectrum end

export injection_shape, minimum_injection_energy
export injection_energy_integral

"""
    injection_shape(spectrum, energy)

Return the dimensionless injection-spectrum shape ``S(E)`` at `energy`.

Concrete spectra must accept an `EnergyUnit` with natural-unit mass dimension
``+1`` and return zero below `minimum_injection_energy(spectrum)`.
"""
function injection_shape end

"""
    minimum_injection_energy(spectrum)

Return the lower boundary of the injection spectrum as an `EnergyUnit` with
natural-unit mass dimension ``+1``.
"""
function minimum_injection_energy end

"""
    injection_energy_integral(spectrum)

Return the energy-normalization integral

```math
\\mathcal{I} = \\int_{E_\\mathrm{min}}^\\infty E\\,S(E)\\,\\mathrm{d}E.
```

The result must be an `EnergyUnit` with natural-unit mass dimension ``+2``.
"""
function injection_energy_integral end
##############################################################################


##############################################################################
export InjectionSpectrum_PowerLawWithExponentialCutoff

"""
    InjectionSpectrum_PowerLawWithExponentialCutoff(
        index,
        minimum_energy,
        cutoff_energy;
        cutoff_index=1,
        reference_energy=EU(1),
    )

Power-law injection shape with a generalized exponential cutoff

```math
S(E) = \\left(\\frac{E}{E_0}\\right)^{-\\gamma}
  \\exp\\left[-\\left(\\frac{E}{E_\\mathrm{cut}}\\right)^\\beta\\right]
  \\Theta(E-E_\\mathrm{min}).
```

The three energies must have natural-unit mass dimension ``+1``. `index` is
``\\gamma`` and the positive `cutoff_index` is ``\\beta``. All energies are
stored in the canonical [`EU`](@ref) basis, and their underlying numeric values
are promoted together with the two indices before conversion to floating
point.
"""
struct InjectionSpectrum_PowerLawWithExponentialCutoff{
    T<:Real,
} <: AbstractInjectionSpectrum
    index::T
    minimum_energy::EU
    cutoff_energy::EU
    cutoff_index::T
    reference_energy::EU

    function InjectionSpectrum_PowerLawWithExponentialCutoff{T}(
        index::T,
        minimum_energy::EU,
        cutoff_energy::EU,
        cutoff_index::T,
        reference_energy::EU,
    ) where {T<:Real}
        _require_finite(index, "injection index")
        _require_dimension(minimum_energy, 1, "minimum injection energy")
        _require_positive_finite(minimum_energy, "minimum injection energy")
        _require_dimension(cutoff_energy, 1, "injection cutoff energy")
        _require_positive_finite(cutoff_energy, "injection cutoff energy")
        _require_positive_finite(cutoff_index, "injection cutoff index")
        _require_dimension(reference_energy, 1, "injection reference energy")
        _require_positive_finite(reference_energy, "injection reference energy")
        return new{T}(
            index,
            minimum_energy,
            cutoff_energy,
            cutoff_index,
            reference_energy,
        )
    end
end

function InjectionSpectrum_PowerLawWithExponentialCutoff(
    index::Real,
    minimum_energy::EnergyUnit,
    cutoff_energy::EnergyUnit;
    cutoff_index::Real=1,
    reference_energy::EnergyUnit=EU(1),
)
    minimum = _canonical_unit(
        minimum_energy,
        1,
        "minimum injection energy",
    )
    cutoff = _canonical_unit(cutoff_energy, 1, "injection cutoff energy")
    reference = _canonical_unit(
        reference_energy,
        1,
        "injection reference energy",
    )
    _require_finite(index, "injection index")
    _require_positive_finite(minimum, "minimum injection energy")
    _require_positive_finite(cutoff, "injection cutoff energy")
    _require_positive_finite(cutoff_index, "injection cutoff index")
    _require_positive_finite(reference, "injection reference energy")

    gamma, E_min, E_cut, beta, E0 = promote(
        index,
        EUval(minimum),
        EUval(cutoff),
        cutoff_index,
        EUval(reference),
    )
    gamma, E_min, E_cut, beta, E0 =
        float(gamma), float(E_min), float(E_cut), float(beta), float(E0)
    return InjectionSpectrum_PowerLawWithExponentialCutoff{typeof(gamma)}(
        gamma,
        EU(E_min),
        EU(E_cut),
        beta,
        EU(E0),
    )
end

minimum_injection_energy(
    spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
) =
    spectrum.minimum_energy

"""
    injection_shape(
        spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
        energy::EnergyUnit,
    )

Evaluate the power law with generalized exponential cutoff

```math
S(E) = \\left(\\frac{E}{E_0}\\right)^{-\\gamma}
  \\exp\\left[-\\left(\\frac{E}{E_\\mathrm{cut}}\\right)^\\beta\\right]
  \\Theta(E-E_\\mathrm{min}).
```

Return zero below `minimum_injection_energy(spectrum)`. The result contains no
source luminosity or total-energy normalization.
"""
function injection_shape(
    spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
    energy::EnergyUnit,
)
    energy = _canonical_unit(energy, 1, "energy")
    _require_positive_finite(energy, "energy")
    numeric_zero = zero(
        spectrum.index * float(EUval(energy)),
    )
    energy < spectrum.minimum_energy && return numeric_zero

    log_energy_ratio = log(energy / spectrum.reference_energy)
    cutoff_ratio = energy / spectrum.cutoff_energy
    cutoff_power = cutoff_ratio^spectrum.cutoff_index
    isinf(cutoff_power) && return numeric_zero

    value = exp(-spectrum.index * log_energy_ratio - cutoff_power)
    return _require_nonnegative_finite(value, "injection shape")
end

"""
    injection_energy_integral(
        spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
    )

Evaluate the energy-normalization integral analytically as

```math
\\mathcal{I}
= \\frac{E_\\mathrm{cut}^2}{\\beta}
  \\left(\\frac{E_\\mathrm{cut}}{E_0}\\right)^{-\\gamma}
  \\Gamma\\!\\left(
    \\frac{2-\\gamma}{\\beta},
    \\left(\\frac{E_\\mathrm{min}}{E_\\mathrm{cut}}\\right)^\\beta
  \\right),
```

where ``\\Gamma(a,x)`` is the upper incomplete Gamma function. The positive
lower boundary keeps this expression finite for every finite `index`.
"""
function injection_energy_integral(
    spectrum::InjectionSpectrum_PowerLawWithExponentialCutoff,
)
    spectral_index = spectrum.index
    cutoff_index = spectrum.cutoff_index
    gamma_parameter = (2 - spectral_index) / cutoff_index
    gamma_argument = (
        spectrum.minimum_energy / spectrum.cutoff_energy
    )^cutoff_index
    prefactor = spectrum.cutoff_energy^2 / cutoff_index *
                (spectrum.cutoff_energy /
                 spectrum.reference_energy)^(-spectral_index)
    value = prefactor * gamma(
        gamma_parameter,
        gamma_argument,
    )
    return _require_positive_finite(value, "injection energy integral")
end
##############################################################################
