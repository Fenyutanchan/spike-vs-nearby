### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 23377bb1-3b92-4509-87fa-d9d554ef17e3
import Pkg; Pkg.activate(@__DIR__)

# ╔═╡ ad46f1b1-b5bc-4d0a-bb89-131400d835c2
begin
    using CairoMakie, LaTeXStrings
    using NaturalUnits
    using spike_vs_nearby
end

# ╔═╡ 706578ae-ad7b-11f1-93fc-e3771aa0c10f
#VSCODE-MARKDOWN
md"# Preliminaries"

# ╔═╡ 1357562e-f43c-4114-8491-5f9daba19ef8
set_theme!(theme_latexfonts())

# ╔═╡ 2b6f1cbf-535c-4863-911e-78acdbeaec2d
function flux_value(flux::EnergyUnit)
    flux_unit = inv(GeV(1.0) * NU.m^2 * NU.s) # GeV⁻¹ m⁻² s⁻¹ sr⁻¹
    flux_value = flux / flux_unit
    @check_EU_dimension flux_value 0
    return flux_value
end

# ╔═╡ 2b391ec3-5ef0-486c-9473-68100ce6a776
#VSCODE-MARKDOWN
md"# Data"

# ╔═╡ abc5bd2b-f5d4-4403-9297-8138aa94151b
begin
    channels = (:AMS02_e⁻, :AMS02_e⁺, :AMS02_e⁺e⁻, :DAMPE_e⁺e⁻)
    data = (
        AMS02_e⁻=read_AMS02_electron_flux(),
        AMS02_e⁺=read_AMS02_positron_flux(),
        AMS02_e⁺e⁻=combine_AMS02_electron_positron_flux(),
        DAMPE_e⁺e⁻=read_DAMPE_electron_positron_flux(),
    )

    md"Data loaded."
end

# ╔═╡ b57a8e4b-ccfa-4038-a71a-d964f9efdd11
#VSCODE-MARKDOWN
md"# Plots"

# ╔═╡ a6321efc-8a63-4fea-8512-0b4f7efacfd4
figure = let
    figure = Figure()

    axis = Axis(figure[1, 1];
        xlabel=L"E_\mathrm{TOA}~\left[\mathrm{GeV}\right]",
        ylabel=
            L"E^3 J~\left[\mathrm{GeV}^2 \, \mathrm{m}^{-2} \,
                \mathrm{s}^{-1} \, \mathrm{sr}^{-1}\right]",
        xscale=log10, yscale=log10,
    )

    plot_settings = (
        AMS02_e⁻ =
            (color=:seagreen4, label=L"AMS-02 $e^-$",
                 marker=:circle),
        AMS02_e⁺ =
            (color=:darkorange2, label=L"AMS-02 $e^+$",
                 marker=:utriangle),
        AMS02_e⁺e⁻ =
            (color=:dodgerblue3, label=L"AMS-02 $e^+ + e^-$",
                marker=:rect),
        DAMPE_e⁺e⁻ =
            (color=:mediumpurple4, label=L"DAMPE $e^+ + e^-$",
                marker=:diamond),
    )

    for channel ∈ channels
        channel_data = getfield(data, channel)
        channel_settings = getfield(plot_settings, channel)
        channel_color = channel_settings.color
        channel_marker = channel_settings.marker
        channel_label = channel_settings.label

        energies_GeV = map(EUval(GeV) ∘ (m->m.energy), channel_data)
        energy_minima_GeV = map(EUval(GeV)∘(m->m.energy_min), channel_data)
        energy_maxima_GeV = map(EUval(GeV)∘(m->m.energy_max), channel_data)
        flux_vals = map(flux_value∘(m->m.flux), channel_data)
        flux_syserr_vals = map(flux_value∘(m->m.systematic_error), channel_data)
        flux_staterr_vals = map(flux_value∘(m->m.statistical_error), channel_data)
        flux_allerr_vals = sqrt.(flux_syserr_vals.^2 + flux_staterr_vals.^2)
        weighted_flux_vals = energies_GeV.^3 .* flux_vals
        weighted_flux_allerr_vals = energies_GeV.^3 .* flux_allerr_vals

        errorbars!(axis, energies_GeV, weighted_flux_vals,
                    energies_GeV .- energy_minima_GeV,
                    energy_maxima_GeV .- energies_GeV;
                    direction=:x, color=(channel_color, 0.55),
                    whiskerwidth=0)

        errorbars!(axis, energies_GeV,
                    weighted_flux_vals, weighted_flux_allerr_vals;
                    color=(channel_color, 0.75), whiskerwidth=4)

        scatter!(axis, energies_GeV, weighted_flux_vals;
                    color=channel_color, marker=channel_marker,
                    markersize=7, label=channel_label
                )

    end
    axislegend(axis; position=:rb)


    figure
end

# ╔═╡ Cell order:
# ╟─706578ae-ad7b-11f1-93fc-e3771aa0c10f
# ╠═23377bb1-3b92-4509-87fa-d9d554ef17e3
# ╠═ad46f1b1-b5bc-4d0a-bb89-131400d835c2
# ╠═1357562e-f43c-4114-8491-5f9daba19ef8
# ╠═2b6f1cbf-535c-4863-911e-78acdbeaec2d
# ╟─2b391ec3-5ef0-486c-9473-68100ce6a776
# ╠═abc5bd2b-f5d4-4403-9297-8138aa94151b
# ╟─b57a8e4b-ccfa-4038-a71a-d964f9efdd11
# ╠═a6321efc-8a63-4fea-8512-0b4f7efacfd4
