# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using FytcPlotRegistries

##############################################################################
export project_directory, src_directory

"""
    project_directory

Absolute path to the root of the `spike_vs_nearby` Julia project.

The path is resolved from the installed module location, so it does not depend
on the process working directory.
"""
const project_directory = pathof(@__MODULE__) |> dirname |> dirname

"""
    src_directory

Absolute path to the `src` directory of the Julia project `spike_vs_nearby`.
"""
const src_directory = joinpath(project_directory, "src")
##############################################################################


##############################################################################
export data_directory, external_data_directory, output_data_directory

"""
    data_directory

Absolute path to the `data` directory of the Julia project `spike_vs_nearby`.
"""
const data_directory = joinpath(project_directory, "data")

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

Absolute path to the `plots` directory of the Julia project `spike_vs_nearby`.
"""
const plot_directory = joinpath(project_directory, "plots")

"""
    plot_registry

Plot registry initialized from `plots/PlotRegistry.toml`.
"""
const plot_registry =
    joinpath(plot_directory, "PlotRegistry.toml") |> PlotRegistry
##############################################################################
