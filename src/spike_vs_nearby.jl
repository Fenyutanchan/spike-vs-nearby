# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module spike_vs_nearby

# Load upstream layers before their dependents. Dependencies may point only to
# an earlier layer. Files within the same layer must be independent and are
# included in alphabetical order.

# layer 1
include("paths.jl")
include("units.jl")

# layer 2
include("utils.jl")

# layer 3
include("flux.jl")
include("injection_spectra.jl")
include("transport.jl")

# layer 4
include("pwn.jl")

end # module spike_vs_nearby
