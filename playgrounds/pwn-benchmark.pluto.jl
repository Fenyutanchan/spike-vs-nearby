### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ f71852e2-c4a8-4bc8-b222-e19ffbd4185d
import Pkg; Pkg.activate(@__DIR__)

# ╔═╡ d95f3659-a459-4e62-a326-15b872280521
begin
    using CairoMakie, LaTeXStrings
    using NaturalUnits, FytcUtilities
    using spike_vs_nearby
end

# ╔═╡ f977be0e-3d5b-4190-8f0c-08979c44371a
#VSCODE-MARKDOWN
md"""# Preliminaries"""

# ╔═╡ 842fb369-a171-4f97-a7b5-4d1cf1fc5d4d
theme_latexfonts() |> set_theme!

# ╔═╡ 282b1987-db36-49e6-95d3-f883c4aa67cf
#VSCODE-MARKDOWN
md"""
# Geminga and Monogem benchmarks

This notebook evaluates the present-day ``e^-+e^+`` fluxes from Geminga
(PSR J0633+1746) and Monogem (PSR J0659+1414/B0656+14).

The source records are the Geminga and Monogem objects exported by
`spike_vs_nearby`, both read from the pinned ATNF catalogue artifact. Geminga
uses the catalogue distance `DIST_A`; Monogem uses the parallax `PX`, with
``d\,[\mathrm{kpc}]=(\texttt{PX}\,[\mathrm{mas}])^{-1}``. The ages and
present spin-down luminosities are derived from `F0` and `F1`. The ages are
characteristic spin-down ages rather than independent measurements of the true
source ages. The luminosity calculation adopts
``I=10^{45}\,\mathrm{g\,cm^2}``, the canonical neutron-star moment of inertia
used by ATNF for its derived `EDOT` parameter; it is not a measured property of
either source.

The remaining inputs are explicit modelling assumptions. Both sources use
pure-dipole spin-down with ``t_\mathrm{sd}=12\,\mathrm{kyr}``, conversion
efficiency ``\eta_e=0.10``, a ``100\,\mathrm{TeV}`` exponential cutoff, and
homogeneous Galactic diffusion. The injection indices are ``2.2`` for
Geminga and ``2.0`` for Monogem, following the one-zone interpretation of the
[HAWC analysis](https://arxiv.org/abs/1711.06223). These calculations are
reproducible one-zone benchmarks, not fits of the observed slow-diffusion
halos.
"""

# ╔═╡ eaa2e08d-d844-4675-9ea0-db860ed88ec0
differential_flux_unit = inv(one(GeV) * NU.m^2 * NU.s) # GeV⁻¹ m⁻² s⁻¹ sr⁻¹

# ╔═╡ c99a243e-1a37-40ab-bcef-d10ae263e0a0
pwn_observations = let
    (
        Geminga=(
            catalogue=Geminga_ATNF,
            distance=parse(Float64, Geminga_ATNF["DIST_A"].value) * kpc,
            distance_parameter="DIST_A",
            age=pulsar_characteristic_age(Geminga_ATNF),
            present_luminosity=pulsar_spin_down_luminosity(Geminga_ATNF),
        ),
        Monogem=(
            catalogue=Monogem_ATNF,
            distance=inv(parse(Float64, Monogem_ATNF["PX"].value)) * kpc,
            distance_parameter="PX",
            age=pulsar_characteristic_age(Monogem_ATNF),
            present_luminosity=pulsar_spin_down_luminosity(Monogem_ATNF),
        ),
    )
end

# ╔═╡ 82e7cb91-a75b-4f5c-a71a-b8dd8285e4d6
pwn_benchmarks = let
    common = (
        efficiency=0.10,
        spin_down_timescale=12 * kyr,
        minimum_energy=EU(1),
        cutoff_energy=EU(1.0e5),
        cutoff_index=1.0,
        reference_energy=EU(1),
    )

    (
        Geminga=(; common..., injection_index=2.2),
        Monogem=(; common..., injection_index=2.0),
    )
end

# ╔═╡ f002b049-9ce8-47b9-8355-682217363698
transport_benchmark = (
    diffusion_normalization=3.0e28 * NU.cm^2 / NU.s,
    diffusion_index=1 / 3,
    loss_normalization=1.0e-16 * one(GeV) / NU.s,
    loss_index=2.0,
    reference_energy=one(GeV),
)

# ╔═╡ 1b8e62bb-c473-486d-80ce-9a548b17da4d
pwn_spectra = let
    function spectrum(benchmark)
        InjectionSpectrum_PowerLawWithExponentialCutoff(
            benchmark.injection_index,
            benchmark.minimum_energy,
            benchmark.cutoff_energy;
            cutoff_index=benchmark.cutoff_index,
            reference_energy=benchmark.reference_energy,
        )
    end

    (
        Geminga=spectrum(pwn_benchmarks.Geminga),
        Monogem=spectrum(pwn_benchmarks.Monogem),
    )
end

# ╔═╡ 4b597610-516a-4ea0-ac4d-cd14dcccfdfe
pwn_transport = TransportModel_Homogeneous(
    DiffusionModel_PowerLaw(
        transport_benchmark.diffusion_normalization,
        transport_benchmark.diffusion_index,
        transport_benchmark.reference_energy,
    ),
    EnergyLossModel_PowerLaw(
        transport_benchmark.loss_normalization,
        transport_benchmark.loss_index,
        transport_benchmark.reference_energy,
    ),
)

# ╔═╡ 0748c17e-6faf-4797-a021-7b1e9e37b11a
pwn_sources = let
    function source(spectrum, benchmark, observation)
        PWNSource(
            spectrum,
            benchmark.efficiency,
            observation.distance,
            observation.age;
            present_luminosity=observation.present_luminosity,
            spin_down_timescale=benchmark.spin_down_timescale,
        )
    end

    (
        Geminga=source(
            pwn_spectra.Geminga,
            pwn_benchmarks.Geminga,
            pwn_observations.Geminga,
        ),
        Monogem=source(
            pwn_spectra.Monogem,
            pwn_benchmarks.Monogem,
            pwn_observations.Monogem,
        ),
    )
end

# ╔═╡ b3ad36c0-2f55-4406-83f2-12324d322b22
#VSCODE-MARKDOWN
md"""
## Present-day flux

The calculation covers ``10\,\mathrm{GeV} \leq E \leq
5\,\mathrm{TeV}``. The individual Geminga and Monogem contributions and
their sum are shown. The plotted quantity is
``E^3\Phi_{e^-+e^+}(E)`` in
``\mathrm{GeV}^2\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}\,\mathrm{sr}^{-1}``.
The steradian is dimensionless in the natural-unit conversion. Each dashed
line marks the energy where the cooling time from infinite source energy
equals the corresponding source age; across it the source-energy ceiling
changes from finite to ``+\infty``.
"""

# ╔═╡ 25aff8c7-ff45-4cc7-afbb-2a145c762ae8
pwn_energy_grid = let
    values_GeV = geomspace(10.0, 5.0e3, 80)
    (values=EU.(values_GeV), values_GeV=values_GeV)
end

# ╔═╡ da36414c-98c4-4db4-80ab-ff3b81326acf
pwn_results = let
    function propagated_result(source)
        densities = map(pwn_energy_grid.values) do energy
            pwn_number_density(
                source,
                pwn_transport,
                energy;
                rtol=1.0e-7,
            )
        end
        fluxes = isotropic_differential_flux.(densities)
        flux_values = [flux / differential_flux_unit for flux in fluxes]
        weighted_flux_values = pwn_energy_grid.values_GeV .^ 3 .* flux_values

        (
            densities=densities,
            fluxes=fluxes,
            flux_values=flux_values,
            weighted_flux_values=weighted_flux_values,
        )
    end

    Geminga = propagated_result(pwn_sources.Geminga)
    Monogem = propagated_result(pwn_sources.Monogem)
    total_fluxes = Geminga.fluxes .+ Monogem.fluxes
    total_flux_values = Geminga.flux_values .+ Monogem.flux_values
    (
        Geminga=Geminga,
        Monogem=Monogem,
        total_fluxes=total_fluxes,
        total_flux_values=total_flux_values,
        total_weighted_flux_values=
            pwn_energy_grid.values_GeV .^ 3 .* total_flux_values,
    )
end

# ╔═╡ 48bea5c8-5f3f-47c4-bbac-5a13f7f6c275
pwn_figure = let
    function cooling_break(observation)
        transport_benchmark.reference_energy^2 /
        transport_benchmark.loss_normalization /
        observation.age
    end
    cooling_breaks_GeV = (
        Geminga=cooling_break(pwn_observations.Geminga) / EU(),
        Monogem=cooling_break(pwn_observations.Monogem) / EU(),
    )
    weighted_flux_extrema = extrema(vcat(
        pwn_results.Geminga.weighted_flux_values,
        pwn_results.Monogem.weighted_flux_values,
        pwn_results.total_weighted_flux_values,
    ))
    minimum_exponent = floor(Int, log10(first(weighted_flux_extrema)))
    maximum_exponent = ceil(Int, log10(last(weighted_flux_extrema)))
    ytick_exponents = minimum_exponent:maximum_exponent
    ytick_values = exp10.(ytick_exponents)
    yticks = (
        ytick_values,
        [latexstring("10^{", exponent, "}") for exponent in ytick_exponents],
    )

    figure = Figure()
    axis = Axis(figure[1, 1];
        xlabel=L"E\;[\mathrm{GeV}]",
        ylabel=L"E^3\,\Phi_{e^-+e^+}\;[\mathrm{GeV}^2\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}\,\mathrm{sr}^{-1}]",
        limits=(nothing, (first(ytick_values), last(ytick_values))),
        xminorgridvisible=true,
        xminorticks=IntervalsBetween(9),
        xscale=log10,
        yminorgridvisible=true,
        yminorticks=IntervalsBetween(9),
        yscale=log10,
        yticks=yticks,
    )
    lines!(axis,
        pwn_energy_grid.values_GeV,
        pwn_results.Geminga.weighted_flux_values;
        color=:royalblue3,
        linewidth=2,
        label="Geminga",
    )
    lines!(axis,
        pwn_energy_grid.values_GeV,
        pwn_results.Monogem.weighted_flux_values;
        color=:darkorange2,
        linewidth=2,
        label="Monogem",
    )
    lines!(axis,
        pwn_energy_grid.values_GeV,
        pwn_results.total_weighted_flux_values;
        color=:black,
        linewidth=2.5,
        label="Geminga + Monogem",
    )
    vlines!(axis, [cooling_breaks_GeV.Geminga];
        color=(:royalblue3, 0.65),
        linestyle=:dash,
    )
    vlines!(axis, [cooling_breaks_GeV.Monogem];
        color=(:darkorange2, 0.65),
        linestyle=:dash,
    )
    annotation_height = first(ytick_values) * 1.15
    text!(axis,
        cooling_breaks_GeV.Geminga,
        annotation_height;
        text=L"\tau_{\max}(E)=t_{\mathrm{age}}~\mathrm{(Geminga)}",
        align=(:left, :bottom),
        color=:royalblue3,
        rotation=pi / 2,
    )
    text!(axis,
        cooling_breaks_GeV.Monogem,
        annotation_height;
        text=L"\tau_{\max}(E)=t_{\mathrm{age}}~\mathrm{(Monogem)}",
        align=(:left, :bottom),
        color=:darkorange2,
        rotation=pi / 2,
    )
    axislegend(axis; backgroundcolor=(:white, 0.750), position=:lt)
    figure
end

# ╔═╡ 090fe98b-58e2-4728-a1a7-b6cbddb39237
#VSCODE-MARKDOWN
md"""
## Numerical checks

For each source, the checks below verify the present-luminosity constructor,
positivity and finiteness of the spectrum, the cooling-limited source-energy
ceiling, and the stability of the propagation integral when `rtol` is
tightened from `1e-5` to `1e-7`.
"""

# ╔═╡ eca0dd92-6a8f-4494-8c29-2d8896c86acb
pwn_checks = let
    function source_check(source, observation, result)
        probe_energies_GeV = [10.0, 100.0, 1.0e3]
        rows = map(probe_energies_GeV) do energy_GeV
            energy = EU(energy_GeV)
            density_coarse = pwn_number_density(
                source,
                pwn_transport,
                energy;
                rtol=1.0e-5,
            )
            density_fine = pwn_number_density(
                source,
                pwn_transport,
                energy;
                rtol=1.0e-7,
            )
            flux_fine = isotropic_differential_flux(density_fine)
            source_maximum = source_energy_ceiling(
                pwn_transport,
                energy,
                observation.age,
            )

            (
                energy_GeV=energy_GeV,
                source_energy_ceiling_GeV=source_maximum / EU(),
                E3_flux=energy_GeV^3 *
                        flux_fine / differential_flux_unit,
                relative_integral_change=
                    abs(density_coarse / density_fine - 1),
            )
        end

        (
            present_luminosity_relative_error=abs(
                pwn_spin_down_luminosity(source, observation.age) /
                observation.present_luminosity - 1,
            ),
            all_fluxes_finite=all(isfinite, result.flux_values),
            all_fluxes_positive=all(>(0), result.flux_values),
            probe_rows=rows,
        )
    end

    (
        Geminga=source_check(
            pwn_sources.Geminga,
            pwn_observations.Geminga,
            pwn_results.Geminga,
        ),
        Monogem=source_check(
            pwn_sources.Monogem,
            pwn_observations.Monogem,
            pwn_results.Monogem,
        ),
        all_total_fluxes_finite=all(isfinite, pwn_results.total_flux_values),
        all_total_fluxes_positive=all(>(0), pwn_results.total_flux_values),
    )
end

# ╔═╡ Cell order:
# ╟─f977be0e-3d5b-4190-8f0c-08979c44371a
# ╠═f71852e2-c4a8-4bc8-b222-e19ffbd4185d
# ╠═d95f3659-a459-4e62-a326-15b872280521
# ╠═842fb369-a171-4f97-a7b5-4d1cf1fc5d4d
# ╟─282b1987-db36-49e6-95d3-f883c4aa67cf
# ╠═eaa2e08d-d844-4675-9ea0-db860ed88ec0
# ╠═c99a243e-1a37-40ab-bcef-d10ae263e0a0
# ╠═82e7cb91-a75b-4f5c-a71a-b8dd8285e4d6
# ╠═f002b049-9ce8-47b9-8355-682217363698
# ╠═1b8e62bb-c473-486d-80ce-9a548b17da4d
# ╠═4b597610-516a-4ea0-ac4d-cd14dcccfdfe
# ╠═0748c17e-6faf-4797-a021-7b1e9e37b11a
# ╟─b3ad36c0-2f55-4406-83f2-12324d322b22
# ╠═25aff8c7-ff45-4cc7-afbb-2a145c762ae8
# ╠═da36414c-98c4-4db4-80ab-ff3b81326acf
# ╠═48bea5c8-5f3f-47c4-bbac-5a13f7f6c275
# ╟─090fe98b-58e2-4728-a1a7-b6cbddb39237
# ╠═eca0dd92-6a8f-4494-8c29-2d8896c86acb
