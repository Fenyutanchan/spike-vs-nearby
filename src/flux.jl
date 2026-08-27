# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

##############################################################################
export isotropic_differential_flux

"""
    isotropic_differential_flux(number_density; beta=1)

Convert the differential number density of an isotropic particle population
into its differential flux,

```math
\\Phi(E) = \\frac{\\beta(E)}{4\\pi} n(E),
```

where `beta` is the dimensionless speed ``\\beta=v/c``. It must satisfy
``0 \\leq \\beta \\leq 1`` and defaults to the ultrarelativistic limit
``\\beta=1``.

`number_density` and the returned flux are `EnergyUnit` values with
natural-unit mass dimension ``+2``. The function preserves the particle
population represented by `number_density`; it does not perform a charge or
species decomposition.
"""
function isotropic_differential_flux(
    number_density::EnergyUnit;
    beta::Real=1,
)
    density = _canonical_unit(
        number_density,
        2,
        "differential number density",
    )
    _require_nonnegative_finite(beta, "particle speed fraction")
    beta <= 1 || throw(DomainError(
        beta,
        "particle speed fraction must not exceed one",
    ))

    density_value, speed_fraction = promote(EUval(density), beta)
    density_value, speed_fraction = float(density_value), float(speed_fraction)
    four_pi = oftype(density_value, 4) * pi
    return EU(speed_fraction * density_value / four_pi, 2)
end
##############################################################################
