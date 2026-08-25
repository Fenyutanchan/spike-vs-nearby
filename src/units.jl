# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using NaturalUnits

##############################################################################
export EU

"""
    EU

Canonical energy-unit type used by `spike_vs_nearby`.

All quantities represented by `NaturalUnits.EnergyUnit` are converted to a
`GeV` basis at the public API boundary while retaining their natural-unit mass
dimension.
"""
const EU = GeV
##############################################################################
