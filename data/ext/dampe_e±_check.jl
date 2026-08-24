# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using DelimitedFiles
using Printf

const EXPECTED_HEADER = [
    "energy_min GeV",
    "energy_max GeV",
    "energy_mean GeV",
    "energy_mean_error GeV",
    "acceptance m^2sr",
    "acceptance_error m^2sr",
    "number_of_events",
    "background_fraction percent",
    "background_fraction_error percent",
    "electron_positron_flux m^-2sr^-1s^-1GeV^-1",
    "electron_positron_flux_error_statistical m^-2sr^-1s^-1GeV^-1",
    "electron_positron_flux_error_systematic m^-2sr^-1s^-1GeV^-1",
]

const SUPERSCRIPT_CHAR = Dict(
    '-' => '⁻',
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
)

const DEFAULT_CSV = joinpath(@__DIR__, "dampe_e±.csv")

"""Read and validate the published 38-row DAMPE Table 1 CSV."""
function load_dampe_table(path::AbstractString)
    raw, header = readdlm(path, ',', header=true)
    names = String.(vec(header))

    names == EXPECTED_HEADER || error("Unexpected CSV header in $(abspath(path))")
    size(raw) == (38, 12) || error(
        "Expected 38 data rows and 12 columns, found $(size(raw, 1)) rows and $(size(raw, 2)) columns",
    )

    data = try
        Float64.(raw)
    catch exception
        error("All DAMPE table cells must be numeric: $(sprint(showerror, exception))")
    end

    all(isfinite, data) || error("DAMPE table contains a non-finite value")
    all(data[:, 1] .< data[:, 3]) || error("Each mean energy must exceed its lower bin edge")
    all(data[:, 3] .< data[:, 2]) || error("Each mean energy must lie below its upper bin edge")
    all(data[1:end-1, 2] .== data[2:end, 1]) || error("Energy bins are not contiguous")
    data[1, 1] == 24.0 || error("Unexpected first lower bin edge")
    data[end, 2] == 4570.9 || error("Unexpected final upper bin edge")
    all(isinteger, data[:, 7]) || error("Event counts must be integers")
    all(data[:, 4:12] .>= 0) || error("Uncertainties and measured quantities must be non-negative")

    return data
end

superscript(n::Integer) = join(SUPERSCRIPT_CHAR[c] for c in string(n))

function format_flux(value::Real, statistical_error::Real, systematic_error::Real)
    value > 0 || error("Flux values must be positive")
    exponent = floor(Int, log10(value))
    scale = 10.0^exponent

    return @sprintf(
        "(%.2f ± %.2f ± %.2f)×10%s",
        value / scale,
        statistical_error / scale,
        systematic_error / scale,
        superscript(exponent),
    )
end

function format_rows(data::AbstractMatrix)
    return [
        [
            @sprintf("%.1f – %.1f", data[row, 1], data[row, 2]),
            @sprintf("%.1f ± %.1f", data[row, 3], data[row, 4]),
            @sprintf("%.3f ± %.3f", data[row, 5], data[row, 6]),
            @sprintf("%d", round(Int, data[row, 7])),
            @sprintf("(%.1f ± %.1f)%%", data[row, 8], data[row, 9]),
            format_flux(data[row, 10], data[row, 11], data[row, 12]),
        ] for row in axes(data, 1)
    ]
end

function pad_display(text::AbstractString, width::Integer, alignment::Symbol)
    padding = width - textwidth(text)
    padding >= 0 || error("Cell width is smaller than its contents")

    if alignment === :left
        return text * repeat(" ", padding)
    elseif alignment === :right
        return repeat(" ", padding) * text
    elseif alignment === :center
        left_padding = padding ÷ 2
        return repeat(" ", left_padding) * text * repeat(" ", padding - left_padding)
    end

    error("Unsupported alignment: $alignment")
end

function horizontal_rule(widths, left::Char, middle::Char, right::Char)
    segments = [repeat("─", width + 2) for width in widths]
    return string(left, join(segments, middle), right)
end

function table_row(cells, widths; alignments=fill(:center, length(widths)))
    formatted = [
        " " * pad_display(string(cell), width, alignment) * " " for
        (cell, width, alignment) in zip(cells, widths, alignments)
    ]
    return "│" * join(formatted, "│") * "│"
end

function print_dampe_table(io::IO, data::AbstractMatrix)
    headers = [
        "Energy range (GeV)",
        "⟨E⟩ (GeV)",
        "Acceptance (m²×sr)",
        "Counts",
        "Bkg. fraction",
        "Φ(e⁺+e⁻) ± σₛₜₐₜ ± σₛᵧₛ",
    ]
    rows = format_rows(data)
    widths = [
        maximum(textwidth(row[column]) for row in [headers, rows...]) for
        column in eachindex(headers)
    ]
    body_width = sum(widths .+ 2) + length(widths) - 1
    title = "Table 1. The CRE flux in units of (m⁻² s⁻¹ sr⁻¹ GeV⁻¹), with 1σ statistical and systematic errors."
    textwidth(title) <= body_width || error("Table title does not fit the computed width")

    println(io, "┌", repeat("─", body_width), "┐")
    println(io, "│", pad_display(title, body_width, :center), "│")
    println(io, horizontal_rule(widths, '├', '┬', '┤'))
    println(io, table_row(headers, widths))
    println(io, horizontal_rule(widths, '├', '┼', '┤'))
    for row in rows
        println(io, table_row(row, widths; alignments=[:center, :center, :center, :right, :center, :right]))
    end
    println(io, horizontal_rule(widths, '└', '┴', '┘'))
end

function main(args=ARGS)
    length(args) <= 1 || error("Usage: julia dampe_e±_check.jl [input.csv]")
    csv_path = isempty(args) ? DEFAULT_CSV : abspath(args[1])
    data = load_dampe_table(csv_path)
    print_dampe_table(stdout, data)
    println(
        "✓ Validated $(size(data, 1)) rows × $(size(data, 2)) columns; ",
        "energy bins are contiguous from 24.0 to 4570.9 GeV.",
    )

    return nothing
end

main()
