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
    notebook="playgrounds/background-prefit.pluto.jl",
)'
```

The notebook can also be opened directly with a Pluto-capable editor.

## Notebooks

| Notebook | Purpose |
| --- | --- |
| [`background-prefit.pluto.jl`](background-prefit.pluto.jl) | Jointly prefit the empirical electron and positron backgrounds to AMS-02 and DAMPE flux data. |
| [`pwn-benchmark.pluto.jl`](pwn-benchmark.pluto.jl) | Compute the present-day $e^-+e^+$ fluxes from Geminga and Monogem under a reproducible one-zone homogeneous-transport benchmark. |

### Background prefit

`background-prefit.pluto.jl` implements the high-energy prefit described in
[`background-prefit.tex`](../../notes/background-prefit.tex), fitting the eight
broken-power-law parameters with Turing maximum likelihood.
All complete bins above `E_min` enter the independent Gaussian likelihood,
with statistical and systematic errors added in quadrature. Solar modulation is
neglected. The single-charge flux comes from `spike_vs_nearby`.

```julia
result = fit_background(E_min=GeV(40), strategy=:AMS_DAMPE);
result.parameter_table
result.covariance
result.figure
plot_background_fit(result; bands=false) # Reuse the stored result, without refitting.
```

The notebook saves the configured fit figure as
`plots/background-prefit-<strategy>-<E_min>GeV.pdf` and records its source notebook
in `plots/PlotRegistry.toml`.

`:AMS_DAMPE` uses AMS positrons and DAMPE total flux. `:AMS` uses the separate
AMS electron and positron spectra. The configured threshold is `GeV(40)`.

The notebook has four parts: physical model and units, data and parameter domains,
maximum likelihood and prediction uncertainties, and plotting. All eight parameters
are free. The distributions in `background_parameter_domains(data)` define the
allowed ranges and Turing's coordinate transformations. Their densities do not
enter maximum likelihood. Normalizations and index changes are positive, and
the spectral indices are real.

Both break-energy lower limits are 100 GeV. The electron upper limit is the last
AMS electron bin edge (1400 GeV) for `:AMS`, or the last DAMPE bin edge
(4570.9 GeV) for `:AMS_DAMPE`. The positron upper limit is always the last AMS
positron bin edge (1000 GeV). `result.domains` records these distributions.

Four explicit starting tuples are tried with Turing's L-BFGS backend and a scaled
initial line-search step. The converged run with the lowest chi-square is retained.
If no run converges, the best finite result is returned with `converged=false`.
`result.runs`, `best_run_index`, `termination_reason`, and the native `raw_fit`
retain the numerical results. The default iteration limit is 3000 and can be set
with `fit_background(E_min=GeV(50), solver_options=(maxiters=5000,))`.

`parameters` contains the best-fit physical values. The curves and residuals use
these values. The result reports total `χ²` and per-channel `channel_χ²`, with
`dof = npoints - nfree` and `reduced_χ² = χ² / dof`.
The free-parameter count `nfree` is derived from the parameter metadata.
Turing's `vcov` supplies the inverse full Hessian of the negative log likelihood
in the original numerical parameters, without rescaling by reduced chi-square.
The returned `covariance` follows `parameter_names`, in products of `parameter_units`.

`Measurements.correlated_values` constructs the eight correlated parameters from
this matrix. Evaluating the existing physical flux functions propagates their
uncertainties, including correlations between the charge contributions.
`parameter_table` shows local standard errors and the plots show pointwise ±1σ
model bands, without added measurement noise. These are local linear error
estimates. `covariance_status` records whether they are usable. An unconverged fit,
an active parameter bound, or an unusable covariance suppresses the errors and bands.
The residual panel shows `(data - fit) / data_error`, with the model error band
centred at zero in the same units. Plotting reuses the stored `spectra` and `residuals`.
On the logarithmic flux panel, band segments with nonpositive lower edges are
omitted. The numerical errors and the residual panel retain the full local errors.

The full-data preview remains in `data-plot.pluto.jl`. The inference stack is
confined to the playground environment. The implementation has no separate
fixed-parameter or generic fitting framework.

The final AMS-only diagnostic at `GeV(40)` compares electron index-change upper
bounds of 8, 16, and 32. It displays the fitted index change, break energy,
chi-square, bound activity, and convergence in a summary table.

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
