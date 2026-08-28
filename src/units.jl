# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using NaturalUnits: GeV, NaturalUnit

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

##############################################################################
export erg, pc, yr
const erg = 1e-7 * NU.J
const pc = 3.08567758149e16 * NU.m
const yr = 31556925.1 * NU.s

export kpc, kyr
const kpc = 1e3 * pc
const kyr = 1e3 * yr
##############################################################################
