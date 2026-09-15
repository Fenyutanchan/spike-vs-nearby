### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 97dd5ef4-6bb5-5531-829d-7352f2044073
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ 7a26bd59-1ccd-4c2b-a1b8-6d85be4d3935
begin
    using CairoMakie, LaTeXStrings
    theme_latexfonts() |> set_theme!

    md"`Makie` initialized."
end

# ╔═╡ c2162faf-355b-57a8-9d79-52e62f03d31c
begin
    using Turing
    using StatsBase: coef, coefnames, vcov

    using Optim: LBFGS, converged, termination_code
    using LineSearches: InitialStatic
    using LinearAlgebra: Diagonal, Symmetric, isposdef, SingularException, PosDefException

    import Measurements

    md"`Turing` initialized."
end

# ╔═╡ 9c418312-9be5-40ba-8cdd-ebf9da8cff9e
begin
    using NaturalUnits
    using FytcPlotRegistries
    using spike_vs_nearby

    md"My Julia packages loaded."
end

# ╔═╡ 03f0c06c-127b-4ea3-8934-a82b0430a838
script_filename = replace(@__FILE__, r"#==#.*$" => "")

# ╔═╡ 1955b0f1-b919-518c-8c5b-43e2f7b5a6dc
#VSCODE-MARKDOWN
md"""
# Background prefit

`fit_background(; E_min, strategy=:AMS_DAMPE)` fits the eight background parameters and returns their best-fit values, local standard errors, covariance, and fit/residual figure.
`:AMS_DAMPE` uses AMS positrons and DAMPE total flux.
`:AMS` uses AMS electrons and positrons.
All complete bins above `E_min` are retained.

Both charges use the package's broken power law with reference energy ``1 \, \mathrm{GeV}``.
Solar modulation is neglected.
Turing provides maximum-likelihood estimation and the full Hessian covariance.
Measurements propagates the correlated parameter errors to the predicted flux.
The full-data preview is in `data-plot.pluto.jl`.
"""

# ╔═╡ 0ed3e31c-7a52-56d2-98b8-ef905086afe3
#VSCODE-MARKDOWN
md"## 1. Physical model and parameter units"

# ╔═╡ 8f0235b8-7f17-5e00-9941-e93fb48933b2
const background_parameter_metadata = (
    C₋ = (unit_str="GeV^-1 m^-2 s^-1 sr^-1", scale=differential_flux_unit),
    Ebr₋ = (unit_str="GeV", scale=GeV(1.0)),
    γ₋ = (unit_str="", scale=1.0),
    Δγ₋ = (unit_str="", scale=1.0),

    C₊ = (unit_str="GeV^-1 m^-2 s^-1 sr^-1", scale=differential_flux_unit),
    Ebr₊ = (unit_str="GeV", scale=GeV(1.0)),
    γ₊ = (unit_str="", scale=1.0),
    Δγ₊ = (unit_str="", scale=1.0),
)

# ╔═╡ ac32361e-2866-5e20-8d43-907b8a049b2e
begin
    # Units are removed only at the fit interface and restored for the flux API.
    function background_numeric_parameters(param)
        return map(param, background_parameter_metadata) do value, metadata
            ratio = value / metadata.scale
            @check_EU_dimension ratio 0
            EUval(ratio)
        end
    end

    background_physical_parameters(param) = map(param, background_parameter_metadata) do value, metadata
        scale = metadata.scale
        value * (scale isa EnergyUnit ? convert(EU, scale) : scale)
    end

    function group_background_parameters(param)
        return (
            electron=(C=param.C₋, Ebr=param.Ebr₋, γ=param.γ₋, Δγ=param.Δγ₋),
            positron=(C=param.C₊, Ebr=param.Ebr₊, γ=param.γ₊, Δγ=param.Δγ₊),
        )
    end
end

# ╔═╡ 3e1bc7d7-2d87-4a18-8d54-a1abdd58b57e
function background_channel_flux(energy::EnergyUnit, parameters, channel)
    channel == :electron && return background_flux(energy, parameters.electron)
    channel == :positron && return background_flux(energy, parameters.positron)
    channel == :electron_positron &&
        return background_flux(energy, parameters.electron) + background_flux(energy, parameters.positron)
    throw(ArgumentError("Unknown flux channel $channel"))
end

# ╔═╡ d36baa99-5336-5c88-b224-471bfd5ea85d
function background_predictions(physical_parameters, data)
    parameters = group_background_parameters(physical_parameters)
    return NamedTuple{keys(data)}(
        map(keys(data)) do channel
            [background_channel_flux(E, parameters, channel) for E ∈ data[channel].energy]
        end
    )
end

# ╔═╡ 92fdfccd-5471-591f-b2db-60788043e605
#VSCODE-MARKDOWN
md"""
## 2. Data, parameter domains and initial values

The likelihood uses the published representative energies and statistical and systematic errors added in quadrature.
The distributions below specify parameter domains in the units listed above.
Their densities do not enter maximum likelihood.
Normalizations and index changes are positive, and the indices are real.
Both break energies start at ``100 \, \mathrm{GeV}``.
Their upper limits follow the selected AMS electron or DAMPE bins for electrons, and the AMS positron bins for positrons.

The same four initial tuples are tried for each fit.
The converged result with the smallest ``\chi^2`` is retained.
L-BFGS uses a scaled initial line-search step.
"""

# ╔═╡ ce6d9c40-aedd-5ad5-8595-818b8db807c5
function select_background_channel(measurements, E_min)
    @check_EU_dimension E_min 1
    isfinite(EUval(E_min)) && E_min > electron_mass ||
        throw(ArgumentError("Require a finite total-energy threshold above m_e"))
    indices = findall(m -> m.energy_min >= E_min, measurements)
    isempty(indices) && throw(ArgumentError("No complete bins above $E_min"))
    selected = measurements[indices]
    return (
        bin_indices = indices,
        energy = [m.energy for m ∈ selected],
        energy_max = [m.energy_max for m ∈ selected],
        flux = [m.flux for m ∈ selected],
        uncertainty = [sqrt(m.statistical_error^2 + m.systematic_error^2) for m ∈ selected],
    )
end

# ╔═╡ c2defc4c-c686-57a7-acba-4c217c367bc5
function load_background_data(strategy, E_min)
    strategies = (AMS=(:electron, :positron), AMS_DAMPE=(:positron, :electron_positron))
    strategy ∈ keys(strategies) || throw(ArgumentError("Use :AMS or :AMS_DAMPE"))
    readers = (
        electron = read_AMS02_electron_flux,
        positron = read_AMS02_positron_flux,
        electron_positron = read_DAMPE_electron_positron_flux
    )
    channels = strategies[strategy]
    return NamedTuple{channels}(map(channel -> select_background_channel(readers[channel](), E_min), channels))
end

# ╔═╡ 58f03616-ac8d-5c17-822c-5556a8f32ca7
function background_parameter_domains(data)
    upper_limit(channel) = maximum(data[channel].energy_max) |> EUval(GeV)
    electron_channel = haskey(data, :electron) ?
        :electron : :electron_positron
    return (
        C₋=LogNormal(),
        Ebr₋=Uniform(100.0, upper_limit(electron_channel)),
        γ₋=Normal(),
        Δγ₋=LogNormal(),

        C₊=LogNormal(),
        Ebr₊=Uniform(100.0, upper_limit(:positron)),
        γ₊=Normal(),
        Δγ₊=LogNormal(),
    )
end

# ╔═╡ 4ff905dc-f2ef-5ec2-ab25-59bb5f483e51
const background_initial_parameters = [
    (
        C₋ = 200.0 * differential_flux_unit,
        Ebr₋ = GeV(150.0),
        γ₋ = 3.0,
        Δγ₋ = 1.0,

        C₊ = 3.0 * differential_flux_unit,
        Ebr₊ = GeV(150.0),
        γ₊ = 2.7,
        Δγ₊ = 1.0,
    ),
    (
        C₋ = 300.0 * differential_flux_unit,
        Ebr₋ = GeV(500.0),
        γ₋ = 3.1,
        Δγ₋ = 1.0,

        C₊ = 5.0 * differential_flux_unit,
        Ebr₊ = GeV(500.0),
        γ₊ = 2.8,
        Δγ₊ = 2.0,
    ),
    (
        C₋ = 500.0 * differential_flux_unit,
        Ebr₋ = GeV(900.0),
        γ₋ = 3.2,
        Δγ₋ = 2.0,

        C₊ = 8.0 * differential_flux_unit,
        Ebr₊ = GeV(900.0),
        γ₊ = 2.9,
        Δγ₊ = 3.0,
    ),
    (
        C₋ = 250.0 * differential_flux_unit,
        Ebr₋ = GeV(800.0),
        γ₋ = 3.0,
        Δγ₋ = 3.0,

        C₊ = 6.0 * differential_flux_unit,
        Ebr₊ = GeV(800.0),
        γ₊ = 2.8,
        Δγ₊ = 1.5,
    ),
]

# ╔═╡ 933450cd-0ee2-5fc2-add2-63ebe9da9efa
const background_fit_solver = LBFGS(alphaguess=InitialStatic(scaled=true))

# ╔═╡ e9deb9f2-00c3-504c-834f-7655fa61be56
#VSCODE-MARKDOWN
md"""
## 3. Maximum likelihood and prediction uncertainties

Turing maximizes the Gaussian data likelihood, equivalent to minimizing ``\chi^2``.
Its `vcov` returns the inverse full Hessian of ``\chi^2/2`` in the original model parameters.
The experimental errors set its scale, without rescaling by reduced chi-square.
Matrix entries use products of the parameter units listed above.

Measurements propagates this full covariance through the physical flux functions.
Spectrum bands show pointwise local ``1\,\sigma`` errors, including charge correlations and without added measurement noise.
Unconverged fits, active parameter bounds, or an unusable covariance leave these errors unavailable.
The selected run's stopping reason and all starting runs remain accessible in the result.
"""

# ╔═╡ 8ea2a78d-c3d2-4433-aec7-4dbf0cd65789
@model function background_model(observed, errors, data, domains)
    C₋ ~ domains.C₋
    Ebr₋ ~ domains.Ebr₋
    γ₋ ~ domains.γ₋
    Δγ₋ ~ domains.Δγ₋

    C₊ ~ domains.C₊
    Ebr₊ ~ domains.Ebr₊
    γ₊ ~ domains.γ₊
    Δγ₊ ~ domains.Δγ₊

    parameters = background_physical_parameters((; C₋, Ebr₋, γ₋, Δγ₋, C₊, Ebr₊, γ₊, Δγ₊))
    predictions = vcat(
        (values ∘ map)(channel -> flux_value.(channel),
            background_predictions(parameters, data)
        )...
    )
    observed ~ MvNormal(predictions, Diagonal(abs2.(errors)))
end

# ╔═╡ 14f2f2d6-59dc-4c0c-90bb-f23b17b5e06b
function background_mode_covariance(mode, numeric, domains, success)
    success || return (matrix=nothing, status=:not_converged)
    native = try
        matrix = Matrix(vcov(mode))
        (matrix + transpose(matrix)) / 2
    catch err
        err isa SingularException || err isa PosDefException || rethrow()
        return (matrix=nothing, status=:unavailable)
    end
    all(isfinite, native) && isposdef(Symmetric(native)) ||
        return (matrix=nothing, status=:unavailable)
    native_names = Symbol.(string.(coefnames(mode)))
    order = [findfirst(==(name), native_names) for name ∈ keys(numeric)]
    at_bound = any(keys(numeric)) do name
        tolerance = 1e-6 * max(abs(numeric[name]), 1.0)
        min(numeric[name] - minimum(domains[name]), maximum(domains[name]) - numeric[name]) <= tolerance
    end
    return (
        matrix = native[order, order],
        status = at_bound ? :at_bound : :ok
    )
end

# ╔═╡ 931453a6-2122-5afc-8b16-04c39b93ce84
function background_spectrum_band(physical_parameters, covariance, energies, channel)
    param = if isnothing(covariance)
        physical_parameters
    else
        numeric = background_numeric_parameters(physical_parameters)
        uncertain = Measurements.correlated_values(collect(values(numeric)), covariance)
        background_physical_parameters(NamedTuple{keys(numeric)}(Tuple(uncertain)))
    end
    parameters = group_background_parameters(param)
    predictions = [flux_value(background_channel_flux(E, parameters, channel)) for E ∈ energies]
    central = Measurements.value.(predictions) .* differential_flux_unit
    isnothing(covariance) && return (central=central, standard_error=nothing, lower=nothing, upper=nothing)
    standard_error = Measurements.uncertainty.(predictions) .* differential_flux_unit
    return (central=central, standard_error=standard_error,
            lower=central .- standard_error, upper=central .+ standard_error)
end

# ╔═╡ 39b7f345-0f20-5692-944c-d491b90c175a
function parameter_table(rows)
    labels = (
        C₋=raw"C_-", Ebr₋=raw"E_{\mathrm{br},-}", γ₋=raw"\gamma_-", Δγ₋=raw"\Delta\gamma_-",
        C₊=raw"C_+", Ebr₊=raw"E_{\mathrm{br},+}", γ₊=raw"\gamma_+", Δγ₊=raw"\Delta\gamma_+",
    )
    number(value) = ismissing(value) ? raw"\text{unavailable}" :
        replace(string(round(value; sigdigits=5)), r"e([+-]?\d+)$" => s"\\times 10^{\1}")
    unit(text) = isempty(text) ? "1" :
        raw"\mathrm{" * replace(text, r"\^(-?\d+)" => s"^{\1}", " " => raw"\,") * "}"
    header = raw"| Parameter | Best fit | Local ``1 \, \sigma`` error | Unit |" * "\n|---|---:|---:|---|\n"
    lines = [
        "| ``$(labels[row.parameter])`` | ``$(number(row.value))`` | ``$(number(row.standard_error))`` | ``$(unit(row.unit_str))`` |" for row ∈ rows
    ]
    return Markdown.parse(header * join(lines, '\n'))
end

# ╔═╡ 2d16ac6f-06fc-5e1a-80ca-e34ddc649a09
#VSCODE-MARKDOWN
md"""
## 4. Fit and residual plots

`plot_background_fit(result; bands=true)` redraws the stored result without refitting.
The upper panel shows ``E^3 J`` and its local ``1 \, \sigma`` model band.
The lower panel shows ``(J_\mathrm{data} - J_\mathrm{fit}) / \sigma`` with unit data error bars.
Its coloured strip is ``\pm\sigma_\mathrm{model}/\sigma``, showing model uncertainty about the fitted curve.
The grey strip marks ``\pm 1``.
The logarithmic flux panel shows the band only where its lower edge is positive.
"""

# ╔═╡ 490fbbe5-ec8d-4d97-a037-09e85221f88a
function plot_background_fit(result; bands=true)
    function log_ticks(low, high)
        exponents = (floor(Int, log10(low)) - 1):(ceil(Int, log10(high)) + 1)
        return (exp10.(exponents), [latexstring("10^{", exponent, "}") for exponent ∈ exponents])
    end
    figure = Figure(size=(1050, 760), fontsize=17)
    show_bands = bands && result.covariance_status == :ok
    channel_labels = (electron=L"AMS $e^-$", positron=L"AMS $e^+$", electron_positron=L"DAMPE $e^-+e^+$")
    channel_summary = join([
        L"%$(channel_labels[channel]): $\chi^2=%$(round(result.channel_χ²[channel]; digits=2)),\ N=%$(length(result.data[channel].flux))$"
        for channel ∈ keys(result.data)
    ], "   |   ")
    strategy_label = replace(string(result.strategy), "_" => " + ")
    fit_summary = L"%$(strategy_label), $E_{\min}=%$(EUval(GeV, result.E_min))\,\mathrm{GeV},\quad \chi^2/\mathrm{dof}=%$(round(result.reduced_χ²; sigdigits=5))$"
    summary_layout = GridLayout()
    figure[0, 1] = summary_layout
    Label(summary_layout[1, 1], fit_summary; fontsize=14, tellwidth=false)
    Label(summary_layout[2, 1], latexstring(channel_summary); fontsize=14, tellwidth=false)
    rowgap!(summary_layout, 2)
    flux_axis = Axis(figure[1, 1]; title="Fit and retained measurements",
        ylabel=L"E^3 J\;[\mathrm{GeV}^2\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}\,\mathrm{sr}^{-1}]",
        xscale=log10, yscale=log10,
        xticks=log_ticks, yticks=log_ticks,
        xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
        xminorticksvisible=true, yminorticksvisible=true)
    pull_axis = Axis(figure[2, 1];
        xlabel=L"E\;[\mathrm{GeV}]",
        ylabel=L"(J_\mathrm{data} - J_\mathrm{fit}) / \sigma",
        xscale=log10,
        xticks=log_ticks, xminorticks=IntervalsBetween(9), xminorticksvisible=true,
        title=show_bands ? L"Grey: $\pm1$. Coloured: $\pm\sigma_{\mathrm{model}}/\sigma$." : L"Grey: $\pm1$.",
        titlesize=13
    )
    colors = (electron=:seagreen4, positron=:darkorange2, electron_positron=:mediumpurple4)
    markers = (electron=:circle, positron=:utriangle, electron_positron=:diamond)
    energies = EUval.(Ref(GeV), result.spectra.energy)
    hspan!(pull_axis, [-1.0], [1.0]; color=(:grey, 0.12))
    for channel ∈ keys(result.data)
        measurement, prediction = result.data[channel], result.spectra[channel]
        residual, color = result.residuals[channel], colors[channel]
        displayed = energies .<= EUval(GeV, maximum(measurement.energy_max))
        grid = energies[displayed]
        data_energies = EUval.(Ref(GeV), measurement.energy)
        if show_bands
            lower, upper = flux_value.(prediction.lower[displayed]), flux_value.(prediction.upper[displayed])
            positive = lower .> 0
            if any(positive)
                # Leave gaps where the local Gaussian band crosses zero on the log axis.
                band!(flux_axis, grid, grid.^3 .* ifelse.(positive, lower, NaN),
                      grid.^3 .* ifelse.(positive, upper, NaN); color=(color, 0.17))
            end
            band!(pull_axis, data_energies, residual.model_lower, residual.model_upper; color=(color, 0.17))
        end
        lines!(flux_axis, grid, grid.^3 .* flux_value.(prediction.central[displayed]);
            color = color,
            linewidth = 2,
            label = channel_labels[channel]
        )
        errorbars!(flux_axis,
            data_energies, data_energies.^3 .* flux_value.(measurement.flux),
            data_energies.^3 .* flux_value.(measurement.uncertainty);
            color = (color, 0.65)
        )
        scatter!(flux_axis, data_energies, data_energies.^3 .* flux_value.(measurement.flux);
            color = color,
            marker = markers[channel],
            markersize = 6
        )
        errorbars!(pull_axis, data_energies, residual.central, ones(length(data_energies));
            color = (color, 0.55),
            whiskerwidth = 3
        )
        scatter!(pull_axis, data_energies, residual.central;
            color = color,
            marker = markers[channel],
            markersize = 6
        )
    end
    hlines!(pull_axis, [0.0]; color=(:black, 0.6), linestyle=:dash)
    xlims!(flux_axis, first(energies), last(energies))
    linkxaxes!(flux_axis, pull_axis)
    hidexdecorations!(flux_axis; grid=false)
    axislegend(flux_axis; position=:lb)
    rowsize!(figure.layout, 1, Relative(0.65))
    return figure
end

# ╔═╡ d6effab7-80e6-57ff-aaee-028806ecdd65
"""
Fit the eight-parameter background and return its covariance and fit/residual figure.
`domain_overrides` replaces named parameter domains for bound comparisons.
"""
function fit_background(; E_min,
    strategy = :AMS_DAMPE,
    domain_overrides::NamedTuple = (;),
    solver_options = (;)
)
    data = load_background_data(strategy, E_min)
    domains = merge(background_parameter_domains(data), domain_overrides)
    observed = vcat((flux_value.(channel.flux) for channel ∈ values(data))...)
    errors = vcat((flux_value.(channel.uncertainty) for channel ∈ values(data))...)
    names = keys(background_parameter_metadata)
    nfree = length(names)
    dof = length(observed) - nfree
    dof > 0 || throw(ArgumentError("The fit needs more than $nfree retained data points"))

    settings = merge((maxiters=3000,), solver_options)
    model = background_model(observed, errors, data, domains)
    runs = map(background_initial_parameters) do initial
        start = background_numeric_parameters(initial)
        try
            mode = maximum_likelihood(model, background_fit_solver;
                initial_params=InitFromParams(start), settings...)
            native_names = Tuple(Symbol.(string.(coefnames(mode))))
            numeric = NamedTuple{names}(NamedTuple{native_names}(Tuple(coef(mode))))
            predictions = background_predictions(background_physical_parameters(numeric), data)
            predicted_flux = vcat((flux_value.(flux) for flux ∈ values(predictions))...)
            (raw_fit=mode, numeric=numeric, χ²=sum(abs2, (observed .- predicted_flux) ./ errors),
             converged=converged(mode.optim_result.original),
             termination_reason=string(termination_code(mode.optim_result.original)), error=nothing)
        catch err
            err isa DomainError || err isa OverflowError || rethrow()
            (raw_fit=nothing, numeric=nothing, χ²=Inf, converged=false,
             termination_reason="EvaluationError", error=sprint(showerror, err))
        end
    end
    finite = findall(run -> isfinite(run.χ²), runs)
    isempty(finite) && error("No finite background fit was obtained: $(first(runs).error)")
    successful = filter(i -> runs[i].converged, finite)
    candidates = isempty(successful) ? finite : successful
    best_index = candidates[argmin([runs[i].χ² for i ∈ candidates])]
    best = runs[best_index]
    covariance = background_mode_covariance(best.raw_fit, best.numeric, domains, best.converged)
    band_covariance = covariance.status == :ok ? covariance.matrix : nothing
    parameters = background_physical_parameters(best.numeric)
    rows = map(enumerate(names)) do (i, name)
        (parameter=name, unit_str=background_parameter_metadata[name].unit_str,
         value=best.numeric[name],
         standard_error=isnothing(band_covariance) ? missing : sqrt(band_covariance[i, i]))
    end

    upper = EUval(GeV, maximum(maximum(channel.energy_max) for channel ∈ values(data)))
    energies = exp.(range(log(EUval(GeV, E_min)), log(upper); length=160))
    energies[1], energies[end] = EUval(GeV, E_min), upper
    energy_grid = GeV.(energies)
    channels = (:electron, :positron, :electron_positron)
    spectra = merge((energy=energy_grid,), NamedTuple{channels}(map(channels) do channel
        background_spectrum_band(parameters, band_covariance, energy_grid, channel)
    end))
    residuals = NamedTuple{keys(data)}(map(keys(data)) do channel
        measurement = data[channel]
        band = background_spectrum_band(parameters, band_covariance, measurement.energy, channel)
        central = flux_value.(band.central)
        sigma = flux_value.(measurement.uncertainty)
        (central=(flux_value.(measurement.flux) .- central) ./ sigma,
         model_lower=isnothing(band.lower) ? nothing : (flux_value.(band.lower) .- central) ./ sigma,
         model_upper=isnothing(band.upper) ? nothing : (flux_value.(band.upper) .- central) ./ sigma)
    end)
    channel_χ² = map(residual -> sum(abs2, residual.central), residuals)
    χ² = sum(channel_χ²)
    result = (
        strategy=strategy, E_min=E_min, data=data, domains=domains, parameters=parameters,
        parameter_names=names, parameter_units=map(entry -> entry.unit_str, background_parameter_metadata),
        parameter_rows=rows, parameter_table=parameter_table(rows),
        covariance=covariance.matrix, covariance_status=covariance.status,
        converged=best.converged, termination_reason=best.termination_reason, raw_fit=best.raw_fit,
        runs=runs, best_run_index=best_index, npoints=length(observed), nfree=nfree, dof=dof,
        χ²=χ², reduced_χ²=χ²/dof, channel_χ²=channel_χ²,
        spectra=spectra, residuals=residuals,
    )
    return merge(result, (figure=plot_background_fit(result),))
end

# ╔═╡ 8e04c129-6164-563c-9409-83d01056c025
# Run configured analyses.
#VSCODE-MARKDOWN
md"## Run: choose the data combination and lower threshold"

# ╔═╡ 603152f8-64ab-5351-8e32-e638194de6b4
background_fit_settings = (
    strategy=:AMS_DAMPE,
    E_min=GeV(40), # Try GeV(50) for comparison.
)

# ╔═╡ d06b0ff8-e6fb-551b-af65-7a1b237ebdd1
background_fit_result = fit_background(; background_fit_settings...);

# ╔═╡ c829080a-365b-563c-8e0a-9023f2dfbe8a
let result = background_fit_result
    Markdown.parse("""
    **$(result.strategy)**: minimum ``\\chi^2 = $(result.χ²)``,
    nominal ``\\mathrm{d.o.f.} = $(result.dof)``, ``\\chi^2 / \\mathrm{d.o.f.} = $(result.reduced_χ²)``.

    - Converged: **$(result.converged)**.
    - Stop: $(result.termination_reason).
    - Covariance: **$(result.covariance_status)**.
    Errors and bands are local ``1 \\, \\sigma`` estimates.
    """)
end

# ╔═╡ 6a7f2c18-ab5d-5567-813a-3e766ade585f
background_fit_result.parameter_table

# ╔═╡ f6e77c35-ca37-5286-9e1b-c7340324c80f
let result = background_fit_result
    threshold_GeV = EUval(GeV, result.E_min)
    plot_file = joinpath(plot_directory,
        "background-prefit-$(result.strategy)-$(threshold_GeV)GeV.pdf")
    save(plot_file, result.figure)
    plot_register!(plot_registry, plot_file, script_filename;
        description="Background prefit with $(result.strategy), E_min=$(threshold_GeV) GeV: " *
                    "E^3 J(E), local 1-sigma model bands and standardized residuals.")
    result.figure
end

# ╔═╡ c51fdffa-c1ab-5007-8514-5c5fa96f5bb3
#VSCODE-MARKDOWN
md"Full-Hessian covariance. Rows and columns follow the parameter table, with entries in the products of the listed units."

# ╔═╡ 5bdbb224-505a-5379-a1df-c70dccef3ebf
background_fit_result.covariance

# ╔═╡ 22f9b156-46ba-4a5b-928a-4a9878e49a12
#VSCODE-MARKDOWN
md"""
## AMS-only diagnostic: the electron index-change bound

At ``E_{\min} = 40 \, \mathrm{GeV}``, the upper bound on ``\Delta{\gamma_-}`` is set to ``8``, ``16``, and ``32`` in turn.
All eight parameters are refitted from the same `background_initial_parameters`, with all other domains unchanged.

The table reports the fitted electron index change, break energy, total ``\chi^2``, bound activity, and convergence for each upper limit.
"""

# ╔═╡ da312964-5dcb-4bf6-b28c-727cb33981e8
;

# ╔═╡ 9b0b0813-e3e3-4c06-b63e-65e76f156156
let
    ams_bound_test_results = map((8.0, 16.0, 32.0)) do upper
        fit = fit_background(;
            E_min = GeV(40),
            strategy = :AMS,
            domain_overrides = (Δγ₋=truncated(LogNormal(); upper=upper),)
        )
        (upper=upper, fit=fit)
    end
    header = raw"| Upper bound on ``\Delta{\gamma_-}`` | Fitted ``\Delta{\gamma_-}`` | " *
        raw"``E_{\mathrm{br},-}~[\mathrm{GeV}]`` | ``\chi^2`` | At upper bound | Converged |" *
        "\n|---:|---:|---:|---:|:---:|:---:|\n"
    lines = map(ams_bound_test_results) do entry
        fit = entry.fit
        param = fit.parameters
        at_upper = isapprox(param.Δγ₋, entry.upper; rtol=1e-6)
        values = (
            entry.upper,
            round(param.Δγ₋; sigdigits=6),
            round(EUval(GeV, param.Ebr₋); sigdigits=6),
            round(fit.χ²; digits=4),
            at_upper,
            fit.converged
        )
        "| " * join(values, " | ") * " |"
    end
    Markdown.parse(header * join(lines, '\n'))
end

# ╔═╡ Cell order:
# ╟─97dd5ef4-6bb5-5531-829d-7352f2044073
# ╟─7a26bd59-1ccd-4c2b-a1b8-6d85be4d3935
# ╟─c2162faf-355b-57a8-9d79-52e62f03d31c
# ╟─9c418312-9be5-40ba-8cdd-ebf9da8cff9e
# ╠═03f0c06c-127b-4ea3-8934-a82b0430a838
# ╟─1955b0f1-b919-518c-8c5b-43e2f7b5a6dc
# ╟─0ed3e31c-7a52-56d2-98b8-ef905086afe3
# ╠═8f0235b8-7f17-5e00-9941-e93fb48933b2
# ╠═ac32361e-2866-5e20-8d43-907b8a049b2e
# ╠═3e1bc7d7-2d87-4a18-8d54-a1abdd58b57e
# ╠═d36baa99-5336-5c88-b224-471bfd5ea85d
# ╟─92fdfccd-5471-591f-b2db-60788043e605
# ╠═ce6d9c40-aedd-5ad5-8595-818b8db807c5
# ╠═c2defc4c-c686-57a7-acba-4c217c367bc5
# ╠═58f03616-ac8d-5c17-822c-5556a8f32ca7
# ╠═4ff905dc-f2ef-5ec2-ab25-59bb5f483e51
# ╠═933450cd-0ee2-5fc2-add2-63ebe9da9efa
# ╟─e9deb9f2-00c3-504c-834f-7655fa61be56
# ╠═8ea2a78d-c3d2-4433-aec7-4dbf0cd65789
# ╠═14f2f2d6-59dc-4c0c-90bb-f23b17b5e06b
# ╠═931453a6-2122-5afc-8b16-04c39b93ce84
# ╠═39b7f345-0f20-5692-944c-d491b90c175a
# ╠═d6effab7-80e6-57ff-aaee-028806ecdd65
# ╟─2d16ac6f-06fc-5e1a-80ca-e34ddc649a09
# ╠═490fbbe5-ec8d-4d97-a037-09e85221f88a
# ╟─8e04c129-6164-563c-9409-83d01056c025
# ╠═603152f8-64ab-5351-8e32-e638194de6b4
# ╠═d06b0ff8-e6fb-551b-af65-7a1b237ebdd1
# ╟─c829080a-365b-563c-8e0a-9023f2dfbe8a
# ╟─6a7f2c18-ab5d-5567-813a-3e766ade585f
# ╠═f6e77c35-ca37-5286-9e1b-c7340324c80f
# ╟─c51fdffa-c1ab-5007-8514-5c5fa96f5bb3
# ╠═5bdbb224-505a-5379-a1df-c70dccef3ebf
# ╟─22f9b156-46ba-4a5b-928a-4a9878e49a12
# ╠═da312964-5dcb-4bf6-b28c-727cb33981e8
# ╟─9b0b0813-e3e3-4c06-b63e-65e76f156156
