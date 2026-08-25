# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using NaturalUnits: EnergyUnit, EUdim, EUval

function unit_isapprox(
    left::EnergyUnit,
    right::EnergyUnit;
    kwargs...,
)
    EUdim(left) == EUdim(right) || return false
    return isapprox(EUval(EU, left), EUval(EU, right); kwargs...)
end
