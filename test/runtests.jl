# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Test
using spike_vs_nearby

# Unit-aware comparisons shared by the physical tests below.
include("test_utils.jl")

# Check physical formulas, normalization, units and observable predictions.
include("background.jl")
include("data.jl")
include("flux.jl")
include("injection_spectra.jl")
include("pwn.jl")
include("transport.jl")
