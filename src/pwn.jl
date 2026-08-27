# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using QuadGK: quadgk

##############################################################################
export PWNSource

"""
    PWNSource

Continuous, point-like PWN source observed at its present age.

The source combines an injection spectrum, a pure-dipole spin-down history,
an effective conversion efficiency, and the distance and time scales needed
for propagation. Its escaping ``e^-+e^+`` injection power is normalized by

```math
\\int_{E_\\mathrm{min}}^\\infty E\\,q_e(E,t)\\,\\mathrm{d}E =
\\eta_e L_\\mathrm{sd}(t).
```

The time-independent integral returned by `injection_energy_integral` is
cached as `energy_normalization`. Use the birth-luminosity constructor when
``L_\\mathrm{sd}^{(0)}`` is known, or the present-luminosity constructor when
the measured ``L_\\mathrm{sd}(t_\\mathrm{age})`` is the available input.
"""
struct PWNSource{
    S<:AbstractInjectionSpectrum,
    T<:Real,
}
    spectrum::S
    initial_luminosity::EU
    spin_down_timescale::EU
    efficiency::T
    distance::EU
    age::EU
    release_delay::EU
    energy_normalization::EU

    function PWNSource{S,T}(
        spectrum::S,
        initial_luminosity::EU,
        spin_down_timescale::EU,
        efficiency::T,
        distance::EU,
        age::EU,
        release_delay::EU,
        energy_normalization::EU,
        ::Val{:canonical},
    ) where {
        S<:AbstractInjectionSpectrum,
        T<:Real,
    }
        _require_dimension(
            initial_luminosity,
            2,
            "initial spin-down luminosity",
        )
        _require_positive_finite(
            initial_luminosity,
            "initial spin-down luminosity",
        )
        _require_dimension(
            spin_down_timescale,
            -1,
            "spin-down timescale",
        )
        _require_positive_finite(
            spin_down_timescale,
            "spin-down timescale",
        )
        _require_nonnegative_finite(efficiency, "PWN efficiency")
        _require_dimension(distance, -1, "PWN distance")
        _require_positive_finite(distance, "PWN distance")
        _require_dimension(age, -1, "PWN age")
        _require_positive_finite(age, "PWN age")
        _require_dimension(release_delay, -1, "PWN release delay")
        _require_nonnegative_finite(release_delay, "PWN release delay")
        release_delay <= age || throw(DomainError(
            release_delay,
            "PWN release delay must not exceed its age",
        ))
        _require_dimension(
            energy_normalization,
            2,
            "injection energy normalization",
        )
        _require_positive_finite(
            energy_normalization,
            "injection energy normalization",
        )
        return new{S,T}(
            spectrum,
            initial_luminosity,
            spin_down_timescale,
            efficiency,
            distance,
            age,
            release_delay,
            energy_normalization,
        )
    end
end

"""
    PWNSource(
        spectrum,
        initial_luminosity,
        spin_down_timescale,
        efficiency,
        distance,
        age;
        release_delay=EU(0, -1),
    )

Construct a PWN source from its birth spin-down luminosity.

`spectrum` is the injection shape ``\\mathcal{S}_\\mathrm{PWN}(E)``,
`initial_luminosity` is ``L_\\mathrm{sd}^{(0)}``, and
`spin_down_timescale` is ``t_\\mathrm{sd}``. The pure-dipole history is

```math
L_\\mathrm{sd}(t) = L_\\mathrm{sd}^{(0)}
\\left(1 + \\frac{t}{t_\\mathrm{sd}}\\right)^{-2}.
```

`efficiency` is the nonnegative effective conversion efficiency ``\\eta_e``.
`distance`, `age`, and `release_delay` are respectively ``d``,
``t_\\mathrm{age}``, and ``t_\\mathrm{rel}``. The release delay is measured
from pulsar birth and must satisfy
``0 \\leq t_\\mathrm{rel} \\leq t_\\mathrm{age}``.

Luminosity has natural-unit mass dimension ``+2``; the spin-down timescale,
distance, age, and release delay have mass dimension ``-1``. Dimensionful
inputs are stored in the canonical [`EU`](@ref) basis, and the injection
energy normalization is computed and cached automatically.
"""
function PWNSource(
    spectrum::S,
    initial_luminosity::EnergyUnit,
    spin_down_timescale::EnergyUnit,
    efficiency::Real,
    distance::EnergyUnit,
    age::EnergyUnit;
    release_delay::EnergyUnit=_zero_unit_like(age, -1),
) where {S<:AbstractInjectionSpectrum}
    luminosity = _canonical_unit(
        initial_luminosity,
        2,
        "initial spin-down luminosity",
    )
    timescale = _canonical_unit(
        spin_down_timescale,
        -1,
        "spin-down timescale",
    )
    eta = float(efficiency)
    source_distance = _canonical_unit(distance, -1, "PWN distance")
    source_age = _canonical_unit(age, -1, "PWN age")
    delay = _canonical_unit(release_delay, -1, "PWN release delay")
    normalization = injection_energy_integral(spectrum)

    L0, t_sd = promote(EUval(luminosity), EUval(timescale))
    L0, t_sd = float(L0), float(t_sd)
    return PWNSource{S,typeof(eta)}(
        spectrum,
        EU(L0, 2),
        EU(t_sd, -1),
        eta,
        source_distance,
        source_age,
        delay,
        normalization,
        Val(:canonical),
    )
end

"""
    PWNSource(
        spectrum,
        efficiency,
        distance,
        age;
        present_luminosity,
        spin_down_timescale,
        release_delay=EU(0, -1),
    )

Construct a PWN source from its measured present-day spin-down luminosity.

`present_luminosity` is ``L_\\mathrm{sd}(t_\\mathrm{age})`` and
`spin_down_timescale` is ``t_\\mathrm{sd}``. The constructor infers the
birth luminosity using

```math
L_\\mathrm{sd}^{(0)} = L_\\mathrm{sd}(t_\\mathrm{age})
\\left(1 + \\frac{t_\\mathrm{age}}{t_\\mathrm{sd}}\\right)^2.
```

`spectrum` is the injection shape ``\\mathcal{S}_\\mathrm{PWN}(E)``,
`efficiency` is ``\\eta_e``, `distance` is ``d``, `age` is
``t_\\mathrm{age}``, and `release_delay` is ``t_\\mathrm{rel}``. The
efficiency must be nonnegative, while the release delay must satisfy
``0 \\leq t_\\mathrm{rel} \\leq t_\\mathrm{age}``.

Luminosity has natural-unit mass dimension ``+2``; the spin-down timescale,
distance, age, and release delay have mass dimension ``-1``. Dimensionful
inputs are stored in the canonical [`EU`](@ref) basis, and the injection
energy normalization is computed and cached automatically.
"""
function PWNSource(
    spectrum::S,
    efficiency::Real,
    distance::EnergyUnit,
    age::EnergyUnit;
    present_luminosity::EnergyUnit,
    spin_down_timescale::EnergyUnit,
    release_delay::EnergyUnit=_zero_unit_like(age, -1),
) where {S<:AbstractInjectionSpectrum}
    luminosity = _canonical_unit(
        present_luminosity,
        2,
        "present spin-down luminosity",
    )
    source_age = _canonical_unit(age, -1, "PWN age")
    timescale = _canonical_unit(
        spin_down_timescale,
        -1,
        "spin-down timescale",
    )
    _require_positive_finite(luminosity, "present spin-down luminosity")
    _require_positive_finite(source_age, "PWN age")
    _require_positive_finite(timescale, "spin-down timescale")

    initial_luminosity = luminosity * (1 + source_age / timescale)^2
    return PWNSource(
        spectrum,
        initial_luminosity,
        timescale,
        efficiency,
        distance,
        source_age;
        release_delay=release_delay,
    )
end
##############################################################################

##############################################################################
export pwn_spin_down_luminosity

"""
    pwn_spin_down_luminosity(source, elapsed_time)

Return the spin-down luminosity of `source` at `elapsed_time` after pulsar
birth.

The elapsed time must have natural-unit mass dimension ``-1`` and be
nonnegative. The result has mass dimension ``+2``.
"""
function pwn_spin_down_luminosity(
    source::PWNSource,
    elapsed_time::EnergyUnit,
)
    time = _canonical_unit(elapsed_time, -1, "elapsed time")
    _require_nonnegative_finite(time, "elapsed time")
    return source.initial_luminosity *
           (1 + time / source.spin_down_timescale)^(-2)
end
##############################################################################

##############################################################################
export pwn_injection_rate

"""
    pwn_injection_rate(source, energy, time_since_birth)

Return the effective escaping ``e^-+e^+`` injection rate per unit energy,

```math
q_e(E,t)
= \\frac{\\eta_e L_\\mathrm{sd}(t)}{\\mathcal{I}}S(E).
```

In natural units this rate has mass dimension zero and is therefore returned
as a `Real`. It is zero before `source.release_delay` and after `source.age`.
"""
function pwn_injection_rate(
    source::PWNSource,
    energy::EnergyUnit,
    time_since_birth::EnergyUnit,
)
    energy = _canonical_unit(energy, 1, "energy")
    _require_positive_finite(energy, "energy")
    time = _canonical_unit(time_since_birth, -1, "time since pulsar birth")
    _require_nonnegative_finite(time, "time since pulsar birth")
    numeric_zero = zero(
        source.efficiency * EUval(source.energy_normalization),
    )
    time < source.release_delay && return numeric_zero
    time > source.age && return numeric_zero

    normalization = source.efficiency *
                    pwn_spin_down_luminosity(source, time) /
                    source.energy_normalization
    return normalization * injection_shape(source.spectrum, energy)
end
##############################################################################

##############################################################################
export pwn_number_density

"""
    pwn_number_density(source, transport, energy; kwargs...)

Return the present-day differential number density of the combined
``e^-+e^+`` population from a continuous point-like PWN source:

```math
n_\\mathrm{PWN}(E,d)
= \\frac{1}{b(E)}
  \\int_E^{E_\\mathrm{s,max}} \\!\\mathrm{d}E_\\mathrm{s}\\,
  K\\!\\left(\\Delta\\widetilde{\\lambda},d\\right)
  q_e\\!\\left(E_\\mathrm{s},t_\\mathrm{age}-\\tau\\right).
```

The finite source age and optional release delay determine
``E_\\mathrm{s,max}``. The result is an `EnergyUnit` with natural-unit mass
dimension ``+2``, corresponding to number per energy per volume. Keyword
arguments are forwarded to `quadgk`.
"""
function pwn_number_density(
    source::PWNSource,
    transport::TransportModel_Homogeneous{D,L},
    energy::EnergyUnit;
    kwargs...,
) where {
    D<:AbstractDiffusionModel,
    L<:EnergyLossModel_PowerLaw,
}
    observed = _canonical_unit(energy, 1, "observed energy")
    _require_positive_finite(observed, "observed energy")
    active_age = source.age - source.release_delay
    iszero(active_age) && return _zero_unit_like(observed, 2)

    source_maximum = source_energy_ceiling(transport, observed, active_age)
    spectrum_minimum = minimum_injection_energy(source.spectrum)
    source_minimum = observed >= spectrum_minimum ? observed : spectrum_minimum
    source_maximum < source_minimum && return _zero_unit_like(observed, 2)
    source_maximum == source_minimum && return _zero_unit_like(observed, 2)

    lower_value, upper_value = promote(
        float(EUval(source_minimum)),
        float(EUval(source_maximum)),
    )
    integral_value, _ = quadgk(lower_value, upper_value; kwargs...) do value
        source_energy = EU(value)
        cooling = cooling_time(transport, source_energy, observed)
        injection_time = source.age - cooling
        rate = pwn_injection_rate(source, source_energy, injection_time)
        kernel = green_kernel(
            diffusion_loss_integral(transport, source_energy, observed),
            source.distance,
        )
        integrand = rate * kernel
        _require_dimension(integrand, 3, "PWN propagation integrand")
        return EUval(EU, integrand)
    end

    energy_integral = EU(integral_value, 4)
    density = energy_integral / energy_loss_rate(transport, observed)
    return _require_nonnegative_finite(density, "PWN number density")
end

##############################################################################
