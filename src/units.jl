# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using NaturalUnits

##############################################################################
export EU, NU

"""
    EU

Canonical energy-unit type used by `spike_vs_nearby`.

All quantities represented by `NaturalUnits.EnergyUnit` are converted to a
`GeV` basis at the public API boundary while retaining their natural-unit mass
dimension.
"""
const EU = GeV

"""
    NU

Natural-unit conversion table associated with [`EU`](@ref).

For example, `NU.cm`, `NU.s`, and `NU.cm^2 / NU.s` represent a length, a time,
and a diffusion coefficient, respectively, in the canonical `GeV` basis.
"""
const NU = NaturalUnit(EU)
##############################################################################
