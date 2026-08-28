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
    notebook="playgrounds/pwn-benchmark.pluto.jl",
)'
```

The notebook can also be opened directly with a Pluto-capable editor.

## Notebooks

| Notebook | Purpose |
| --- | --- |
| [`pwn-benchmark.pluto.jl`](pwn-benchmark.pluto.jl) | Compute the present-day $e^-+e^+$ fluxes from Geminga and Monogem under a reproducible one-zone homogeneous-transport benchmark. |

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
