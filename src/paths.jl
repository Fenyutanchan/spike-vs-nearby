# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using FytcPlotRegistries

##############################################################################
export code_directory

"""
    code_directory

Absolute path to the root of the `spike_vs_nearby` Julia project.

The path is resolved from the installed module location, so it does not depend
on the process working directory.
"""
const code_directory = pathof(@__MODULE__) |> dirname |> dirname
##############################################################################


##############################################################################
export data_directory, external_data_directory, output_data_directory

"""
    data_directory

Absolute path to the project's `data` directory.
"""
const data_directory = joinpath(code_directory, "data")

"""
    external_data_directory

Absolute path to `data/ext`, which contains immutable external data inputs.
"""
const external_data_directory = joinpath(data_directory, "ext")

"""
    output_data_directory

Absolute path to `data/out`, reserved for generated numerical outputs.

This constant identifies the intended location but does not create the
directory.
"""
const output_data_directory = joinpath(data_directory, "out")
##############################################################################


##############################################################################
export plot_directory, plot_registry

"""
    plot_directory

Absolute path to the project's `plots` directory.
"""
const plot_directory = joinpath(code_directory, "plots")

"""
    plot_registry

Plot registry initialized from `plots/PlotRegistry.toml`.
"""
const plot_registry =
    joinpath(plot_directory, "PlotRegistry.toml") |> PlotRegistry
##############################################################################

