# Playgrounds

This directory contains runnable notebooks and scripts for exploring the
scientific calculations provided by `spike_vs_nearby`.

> [!WARNING]
> For technical reasons, the Pluto notebooks are not currently reproducible
> as independent, self-contained files. They require the complete repository
> checkout and the
> [FytcJuliaRegistry.jl](https://github.com/Fenyutanchan/FytcJuliaRegistry.jl.git)
> Julia registry.

## Environment

The notebooks share the environment declared in
[`Project.toml`](Project.toml). From the package root, instantiate it with:

```sh
julia --project=playgrounds -e 'using Pkg; Pkg.instantiate()'
```

Each Pluto notebook activates this environment automatically and uses the
local `spike_vs_nearby` source tree.

Launch a notebook from an environment in which Pluto is available:

```sh
julia -e 'using Pluto; Pluto.run(
    notebook="playgrounds/background-fit.pluto.jl",
)'
```

The notebook can also be opened directly with a Pluto-capable editor.

## Notebooks

| Notebook | Purpose |
| --- | --- |
| [`background-fit.pluto.jl`](background-fit.pluto.jl) | Jointly prefit the empirical electron and positron backgrounds to AMS-02 and DAMPE flux data. |
| [`pwn-benchmark.pluto.jl`](pwn-benchmark.pluto.jl) | Compute the present-day $e^-+e^+$ fluxes from Geminga and Monogem under a reproducible one-zone homogeneous-transport benchmark. |

### Background prefit

`background-fit.pluto.jl` implements the high-energy prefit in
[`background-prefit.tex`](../../notes/background-prefit.tex).
`fit_background(; strategy, E_min)` uses all complete bins above the supplied
lower threshold, through each channel's last published bin. The two strategies
are `:AMS` (AMS electrons and positrons) and `:AMS_DAMPE` (AMS positrons and the
DAMPE total flux). Solar modulation is neglected in the current calculations,
with its impact left for a later estimate.

The single-charge background spectrum and experimental flux-unit conversion
are provided by `spike_vs_nearby`. The notebook assembles the observed charge
channels with `background_channel_flux` and retains the fitting workflow for review:

| Notebook section | Responsibility |
| --- | --- |
| 1. Units, background spectra and observables | Imported flux API, normalization reference and observed-channel assembly. |
| 2. Fit settings | Eight parameter names, units, numerical bounds and initial values. |
| 3. Parameter coordinates | Conversions and explicitly fixed parameters. |
| 4. Measurements, bin selection and chi-square | Data loading, lower-threshold selection and residuals. |
| 5. Fit at a fixed lower threshold | Optimization from the common starting points. |
| 6. Parameter covariance and spectrum uncertainty | Joint covariance and local spectrum error bands. |
| 7. Parameter table, fitted spectra and residuals | Plot an existing result or supplied parameters. |
| 8. One-call interface | `fit_background` and its returned results. |
| 9. Run a specified lower threshold | Configure the strategy and `E_min`, then display the result. |

The configured threshold is `GeV(40)`, with `GeV(50)` available for comparison.
Each charge has four broken-power-law parameters. The break energies have a
100 GeV lower search bound and no upper bound. All parameters carry their
physical units, and the optimizer uses dimensionless coordinates.

```julia
AMS_fit = fit_background(strategy=:AMS, E_min=GeV(50));
DAMPE_fit = fit_background(strategy=:AMS_DAMPE, E_min=GeV(50));
plot_background_fit(AMS_fit)
plot_background_fit(DAMPE_fit)
```

The results include chi-square summaries, fitted parameters, a joint covariance,
and a figure with spectra and standardized residuals. The optional spectrum
bands are pointwise local 1σ bands. Covariance rows and columns follow the
returned parameter order and units. A rank-deficient fit has no full covariance
estimate. The full-data comparison is in `data-plot.pluto.jl`.

The tests load definitions up to the `# Run configured analyses.` marker without
running the editable example. From the package root:

```sh
julia --startup-file=no --project=playgrounds playgrounds/test_background_prefitting.jl
```

### PWN benchmark

`pwn-benchmark.pluto.jl`:

- reads Geminga and Monogem from the pinned ATNF pulsar catalogue artifact;
- derives their characteristic ages and present spin-down luminosities;
- constructs exponentially cut off power-law injection spectra;
- propagates continuous PWN injection through homogeneous diffusion and
  continuous energy losses;
- evaluates the flux between $10\,\mathrm{GeV}$ and $5\,\mathrm{TeV}$;
- plots the individual and combined $E^3\Phi_{e^-+e^+}(E)$ spectra; and
- reports luminosity reconstruction and integration-convergence diagnostics.

The calculation is a source-only benchmark. Its injection and transport
parameters are explicit modelling assumptions, not a fit to AMS-02, DAMPE, or
the observed Geminga and Monogem slow-diffusion halos.

To use it:

1. Run all cells once to load the ATNF artifact and evaluate the baseline.
2. Edit `pwn_benchmarks` to vary the injection indices, efficiency,
   spin-down timescale, minimum energy, or cutoff parameters.
3. Edit `transport_benchmark` to vary the homogeneous diffusion and
   energy-loss laws.
4. Edit `pwn_energy_grid` to change the plotted energy interval or number
   of sampling points.
5. Inspect `pwn_figure` for the source spectra, `pwn_results` for the
   numerical arrays, and `pwn_checks` for convergence and consistency
   diagnostics.
