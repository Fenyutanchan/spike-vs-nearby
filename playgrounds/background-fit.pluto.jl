### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 97dd5ef4-6bb5-5531-829d-7352f2044073
import Pkg; Pkg.activate(@__DIR__)

# ╔═╡ c2162faf-355b-57a8-9d79-52e62f03d31c
begin
    using NaturalUnits, ForwardDiff, LinearAlgebra, Optim
    using CairoMakie, LaTeXStrings
    using spike_vs_nearby
    set_theme!(theme_latexfonts())
end

# ╔═╡ 1955b0f1-b919-518c-8c5b-43e2f7b5a6dc
md"""
# Background prefit above a lower energy threshold

`fit_background(; strategy, E_min)` fits all available complete bins above `E_min`.
The strategies are `:AMS` (AMS electrons and positrons) and `:AMS_DAMPE`
(AMS positrons and the DAMPE total flux).
Each charge has its own broken power law. Solar modulation is neglected in the
current project calculations. Its impact will be estimated separately later.
All definitions remain in this notebook. The final cells display the fit,
parameter table, covariance, and spectra with standardized residuals.
"""

# ╔═╡ 0ed3e31c-7a52-56d2-98b8-ef905086afe3
md"""
## 1. Units, background spectra and observables

The empirical spectra are compared directly with the measured total-energy
fluxes. DAMPE measures the sum of the electron and positron contributions.
"""

# ╔═╡ 41d8f99a-f021-5112-9158-4036bc61660a
begin
    const background_reference_energy = GeV(1.0)
    const differential_flux_unit = inv(GeV(1.0) * NU.m^2 * NU.s)
    const background_fit_strategies = (AMS=(:electron, :positron),
        AMS_DAMPE=(:positron, :electron_positron))
    const background_channels = (:electron, :positron, :electron_positron)
end

# ╔═╡ 3c6ab0db-9a1c-558d-8338-c1eeb71b1d0b
function flux_value(flux)
    ratio = flux / differential_flux_unit
    @check_EU_dimension ratio 0
    return EUval(ratio)
end

# ╔═╡ b5ce58d3-1cd3-53dd-9a48-937949dea816
begin
    function background_flux(energy::EnergyUnit, param;
                             reference=background_reference_energy)
        @check_EU_dimension energy 1
        @check_EU_dimension reference 1
        @check_EU_dimension param.Ebr 1
        @check_EU_dimension param.C 2
        @check_EU_dimension param.γ 0
        @check_EU_dimension param.Δγ 0
        isfinite(EUval(energy)) && energy > electron_mass ||
            throw(DomainError(energy, "Total energy must exceed the rest energy"))
        param.C > zero(param.C) && param.Ebr > electron_mass && param.Δγ > 0 ||
            throw(ArgumentError("Invalid background parameters"))
        reference > zero(reference) || throw(DomainError(reference))
        # Stable log(1 + exp(x)), including AD through the spectral break.
        x = param.Δγ * log(EUval(energy / param.Ebr))
        softplus = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
        return param.C * exp(-param.γ * log(EUval(energy / reference)) - softplus)
    end

    function background_flux(energy::EnergyUnit, parameters, channel)
        channel == :electron && return background_flux(energy, parameters.electron)
        channel == :positron && return background_flux(energy, parameters.positron)
        channel == :electron_positron &&
            return background_flux(energy, parameters.electron) +
                   background_flux(energy, parameters.positron)
        throw(ArgumentError("Unknown flux channel $channel"))
    end
end

# ╔═╡ 92fdfccd-5471-591f-b2db-60788043e605
md"""
## 2. Fit settings
"""

# ╔═╡ 3744a5cc-61e6-4cc5-a2a9-1bac68a9cb63
md"""
Each named parameter stores its displayed unit, numerical scale and coordinate
choice. Bounds and initial values are matched by name. All eight parameters are
free by default. `fixed` restores any explicitly fixed values in the model.
The numerical settings and starting points can be overridden in each call.
Both break energies have a 100 GeV lower bound and no upper bound.
"""

# ╔═╡ 8f0235b8-7f17-5e00-9941-e93fb48933b2
const background_parameter_metadata = (
    C₋ = (unit_str="GeV^-1 m^-2 s^-1 sr^-1",
            scale=differential_flux_unit, log_coordinate=true),
    Ebr₋ = (unit_str="GeV", scale=GeV(1.0), log_coordinate=true),
    γ₋ = (unit_str="", scale=1.0, log_coordinate=false),
    Δγ₋ = (unit_str="", scale=1.0, log_coordinate=true),

    C₊ = (unit_str="GeV^-1 m^-2 s^-1 sr^-1",
            scale=differential_flux_unit, log_coordinate=true),
    Ebr₊ = (unit_str="GeV", scale=GeV(1.0), log_coordinate=true),
    γ₊ = (unit_str="", scale=1.0, log_coordinate=false),
    Δγ₊ = (unit_str="", scale=1.0, log_coordinate=true),
)

# ╔═╡ 58f03616-ac8d-5c17-822c-5556a8f32ca7
const background_parameter_bounds = (
    lower=(
        C₋=1.0e-4 * differential_flux_unit, Ebr₋=GeV(100.0),
        γ₋=-3.0, Δγ₋=0.001,

        C₊=1.0e-4 * differential_flux_unit, Ebr₊=GeV(100.0),
        γ₊=-3.0, Δγ₊=0.001,
    ),
    upper=(
        C₋=1.0e5 * differential_flux_unit, Ebr₋=GeV(Inf),
        γ₋=6.0, Δγ₋=8.0,

        C₊=1.0e5 * differential_flux_unit, Ebr₊=GeV(Inf),
        γ₊=6.0, Δγ₊=8.0,
    ),
)

# ╔═╡ ff397277-672e-5813-8fe0-a4d6ebc1606e
# Example entry: γ₋=(value=2.0, reason="independent constraint").
const background_fixed_parameters = (;)

# ╔═╡ 4ff905dc-f2ef-5ec2-ab25-59bb5f483e51
const background_initial_parameters = [
    (
        C₋=200.0 * differential_flux_unit, Ebr₋=GeV(150.0),
        γ₋=3.0, Δγ₋=1.0,

        C₊=3.0 * differential_flux_unit, Ebr₊=GeV(150.0),
        γ₊=2.7, Δγ₊=1.0,
    ),
    (
        C₋=300.0 * differential_flux_unit, Ebr₋=GeV(500.0),
        γ₋=3.1, Δγ₋=1.0,

        C₊=5.0 * differential_flux_unit, Ebr₊=GeV(500.0),
        γ₊=2.8, Δγ₊=2.0,
    ),
    (
        C₋=500.0 * differential_flux_unit, Ebr₋=GeV(1500.0),
        γ₋=3.2, Δγ₋=2.0,

        C₊=8.0 * differential_flux_unit, Ebr₊=GeV(1500.0),
        γ₊=2.9, Δγ₊=3.0,
    ),
]

# ╔═╡ 933450cd-0ee2-5fc2-add2-63ebe9da9efa
default_numerics() = (
    iterations=3_000, gradient_tolerance=1.0e-6,
    objective_relative_tolerance=1.0e-12,
    initial_step=0.1, svd_relative_tolerance=1.0e-8,
    bound_tolerance=1.0e-6, max_condition_number=1.0e6,
)

# ╔═╡ 3dd95f9b-4f2c-56cf-9adf-735bee6ea300
md"""
## 3. Parameter coordinates

Physical parameters carry units. The optimizer uses dimensionless coordinates,
with logarithms for the positive normalization, break energy and index change.
Fixed parameters are restored when predictions are evaluated. For example,
`fixed=(γ₋=(value=2.0, reason="independent constraint"),)` fixes one slope.
Methods of each conversion function share a cell so that Pluto can update them
together.
"""

# ╔═╡ 598bba1a-953b-5085-93aa-042c4366f6fb
function background_parameter_value(value, name)
    scale = background_parameter_metadata[name].scale
    @check_EU_dimension value EUdim(scale)
    return EUval(value / scale)
end

# ╔═╡ 01cd9677-d31a-5cc9-849f-f89983ff1668
begin
    background_parameter_values(physical::NamedTuple) = collect(promote((
        background_parameter_value(physical[name], name)
        for name ∈ keys(background_parameter_metadata)
    )...))

    function background_parameter_values(coordinates, space)
        length(coordinates) == length(space.free_indices) ||
            throw(ArgumentError("Coordinate and free-parameter counts differ"))
        physical = convert.(eltype(coordinates), space.template)
        for (j, i) ∈ enumerate(space.free_indices)
            physical[i] = space.metadata[i].log_coordinate ? exp(coordinates[j]) : coordinates[j]
        end
        return physical
    end
end

# ╔═╡ 67e02064-abbe-51c7-bbc9-af8ee56b45d7
begin
    physical_background_parameters(values) = NamedTuple{keys(background_parameter_metadata)}(Tuple(
        values[i] * background_parameter_metadata[name].scale
        for (i, name) ∈ enumerate(keys(background_parameter_metadata))))

    physical_background_parameters(coordinates, space) =
        physical_background_parameters(background_parameter_values(coordinates, space))
end

# ╔═╡ d15f2c11-32bb-59c0-908c-35d9ee1e3a44
function encode_background_parameters(physical, space)
    length(physical) == length(space.names) ||
        throw(ArgumentError("Expected $(length(space.names)) physical parameters"))
    values = background_parameter_values(physical)
    all(isfinite, values) &&
    all(values .>= space.lower .- 8eps(Float64) .* (1 .+ abs.(space.lower))) &&
    all(values .<= space.upper .+ 8eps(Float64) .* (1 .+ abs.(space.upper))) ||
        throw(ArgumentError("Starting parameters lie outside their bounds"))
    bounded = clamp.(values, space.lower, space.upper)
    isempty(space.free_indices) && return Float64[]
    return [
        space.metadata[i].log_coordinate ? log(bounded[i]) : bounded[i]
        for i ∈ space.free_indices
    ]
end

# ╔═╡ ac32361e-2866-5e20-8d43-907b8a049b2e
function unpack_background_parameters(param::NamedTuple)
    return (
        electron=(C=param.C₋, Ebr=param.Ebr₋, γ=param.γ₋, Δγ=param.Δγ₋),
        positron=(C=param.C₊, Ebr=param.Ebr₊, γ=param.γ₊, Δγ=param.Δγ₊),
    )
end

# ╔═╡ fcd83373-6b69-576b-8eff-f440fd700da9
function background_parameter_space(bounds, fixed=(;))
    names = keys(background_parameter_metadata)
    lower, upper = background_parameter_values(bounds.lower), background_parameter_values(bounds.upper)
    length(lower) == length(upper) == length(names) ||
        throw(ArgumentError("Expected $(length(names)) parameter bounds"))
    all(isfinite, lower) && all(x -> isfinite(x) || x == Inf, upper) &&
    all(lower .< upper) ||
        throw(ArgumentError("Bounds must be ordered with finite lower bounds"))
    for param ∈ values(unpack_background_parameters(bounds.lower))
        param.C > zero(param.C) && param.Ebr > electron_mass &&
        param.Δγ > 0 ||
            throw(ArgumentError("Bounds violate the physical parameter domains"))
    end
    all(name -> name ∈ names, keys(fixed)) ||
        throw(ArgumentError("Unknown fixed parameter"))
    # Free entries are replaced by optimizer coordinates. Only fixed entries
    # need stored values, so a finite lower-bound template also supports Inf.
    template = copy(lower)
    for (name, specification) ∈ pairs(fixed)
        index = findfirst(==(name), names)
        value = background_parameter_value(specification.value, name)
        isfinite(value) && lower[index] <= value <= upper[index] ||
            throw(ArgumentError("Fixed value for $name lies outside its bounds"))
        isempty(strip(specification.reason)) &&
            throw(ArgumentError("Record a reason for fixing $name"))
        template[index] = value
    end
    free_indices = findall(name -> !haskey(fixed, name), collect(names))
    encode(value, index) = background_parameter_metadata[index].log_coordinate ? log(value) : value
    return (
        names=names, metadata=background_parameter_metadata,
        lower=lower, upper=upper, fixed=fixed, template=template,
        free_indices=free_indices,
        coordinate_lower=[encode(lower[i], i) for i ∈ free_indices],
        coordinate_upper=[encode(upper[i], i) for i ∈ free_indices],
    )
end

# ╔═╡ 9d5b5516-2a36-5f28-a6a6-12d4dee1e6f3
decode_background_parameters(coordinates, space) =
    unpack_background_parameters(physical_background_parameters(coordinates, space))

# ╔═╡ fef949b0-4cad-54b6-8dd0-4354032c2919
md"""
## 4. Measurements, bin selection and ``\chi^2``

All complete bins with a lower edge at or above `E_min` enter the fit.
Each selected channel contributes through its last published bin.
The model is evaluated at the published representative energies. Statistical and
systematic errors are added in quadrature. The objective sums squared
standardized residuals over the two channels of the chosen strategy.
"""

# ╔═╡ 1f3cbf9f-7420-5569-a632-2aaf97ec9b19
function background_flux_data(measurements::AbstractVector{<:DifferentialFluxMeasurement})
    statistical_error = [convert(GeV, m.statistical_error) for m ∈ measurements]
    systematic_error = [convert(GeV, m.systematic_error) for m ∈ measurements]
    return (
        bin_indices=collect(eachindex(measurements)),
        energy_min=[convert(GeV, m.energy_min) for m ∈ measurements],
        energy_max=[convert(GeV, m.energy_max) for m ∈ measurements],
        energy=[convert(GeV, m.energy) for m ∈ measurements],
        flux=[convert(GeV, m.flux) for m ∈ measurements],
        statistical_error=statistical_error,
        systematic_error=systematic_error,
        uncertainty=sqrt.(abs2.(statistical_error) .+ abs2.(systematic_error)),
    )
end

# ╔═╡ c2defc4c-c686-57a7-acba-4c217c367bc5
"""Read only the two flux channels used by the requested strategy."""
function load_background_data(strategy=:AMS)
    channels = background_fit_strategies[strategy]
    readers = (electron=read_AMS02_electron_flux,
               positron=read_AMS02_positron_flux,
               electron_positron=read_DAMPE_electron_positron_flux)
    return NamedTuple{channels}(Tuple(background_flux_data(readers[channel]())
                                     for channel ∈ channels))
end

# ╔═╡ aeb77dc9-f168-5529-b9ee-eaaced92bf64
function validate_background_threshold(E_min)
    @check_EU_dimension E_min 1
    isfinite(EUval(E_min)) && electron_mass <= E_min ||
        throw(ArgumentError("Require a finite E_min >= m_e"))
    return E_min
end

# ╔═╡ 484285e8-90e0-526b-94ef-47d81e21c674
function background_energy_grid(data, E_min; points=160)
    upper_edges = vcat((channel.energy_max for channel ∈ values(data))...)
    isempty(upper_edges) && throw(ArgumentError("No retained data to plot"))
    lower = max(EUval(GeV, E_min), nextfloat(EUval(GeV, electron_mass)))
    upper = EUval(GeV, maximum(upper_edges))
    return GeV.(exp.(range(log(lower), log(upper); length=points)))
end

# ╔═╡ ce6d9c40-aedd-5ad5-8595-818b8db807c5
function select_background_channel(data, E_min)
    validate_background_threshold(E_min)
    keep = findall(data.energy_min .>= E_min)
    selected = map(values -> values[keep], data)
    all(x -> isfinite(EUval(x)), selected.flux) &&
    all(x -> isfinite(EUval(x)) && x > zero(x), selected.uncertainty) ||
        throw(ArgumentError("Retained fluxes must be finite and errors positive"))
    return selected
end

# ╔═╡ a936d1ed-3d0e-5a53-a934-1c99fe6a469b
select_background_data(data, E_min) =
    map(channel -> select_background_channel(channel, E_min), data)

# ╔═╡ d36baa99-5336-5c88-b224-471bfd5ea85d
function background_predictions(coordinates, data, space)
    parameters = decode_background_parameters(coordinates, space)
    return NamedTuple{keys(data)}(Tuple(
        GeV[convert(GeV, background_flux(E, parameters, channel))
            for E ∈ data[channel].energy]
        for channel ∈ keys(data)
    ))
end

# ╔═╡ 2fc9db5e-1b73-581e-8905-5a1a09be12f4
function background_pulls(coordinates, data, space)
    predictions = background_predictions(coordinates, data, space)
    return map(data, predictions) do channel, prediction
        pulls = Vector{eltype(coordinates)}(undef, length(channel.flux))
        for i ∈ eachindex(channel.flux)
            pull = (channel.flux[i] - prediction[i]) / channel.uncertainty[i]
            @check_EU_dimension pull 0
            pulls[i] = EUval(pull)
        end
        pulls
    end
end

# ╔═╡ 57e0e552-d709-5a89-9939-ebddd38a46c0
background_standardized_residuals(coordinates, data, space) =
    vcat(values(background_pulls(coordinates, data, space))...)

# ╔═╡ c271919b-f8a3-5961-aef4-8fc8f854c0d8
background_chi2(coordinates, data, space) =
    sum(abs2, background_standardized_residuals(coordinates, data, space))

# ╔═╡ e9deb9f2-00c3-504c-834f-7655fa61be56
md"""
## 5. Fit at a fixed lower threshold

Each supplied starting point is optimized within the parameter bounds.
The converged fit with the smallest chi-square is used. If no run converges,
the best finite result is returned with `converged=false`.
The nominal degrees of freedom are the selected-point count minus the number
of free parameters and must be positive.
"""

# ╔═╡ aeae588a-599b-549a-9deb-39b7fd147d4b
function fit_background_from(initial_parameters, data, space, numerics)
    initial = encode_background_parameters(initial_parameters, space)
    objective = z -> background_chi2(z, data, space)
    if isempty(initial)
        return (initial_coordinates=initial, coordinates=initial, chi2=objective(initial),
                converged=true, iterations=0, projected_gradient=0.0)
    end
    # Roundoff in a line search can move a boundary value slightly outside
    # the box. Evaluate at the projected point and retain its model gradient.
    project = z -> clamp.(z, space.coordinate_lower, space.coordinate_upper)
    # A constant numerical rescaling prevents the first quasi-Newton step from
    # jumping into a flat, nearly zero-flux corner. It leaves the minimum and
    # the reported chi² unchanged and does not alter the supplied initial point.
    scale = max(objective(initial), 1.0)
    scaled_objective = z -> objective(z) / scale
    bounded_objective = z -> scaled_objective(project(z))
    gradient! = (g, z) -> ForwardDiff.gradient!(g, scaled_objective, project(z))
    result = Optim.optimize(
        bounded_objective, gradient!, space.coordinate_lower, space.coordinate_upper,
        initial, Optim.LBFGSB(m=15, stepsize=numerics.initial_step),
        Optim.Options(
            iterations=numerics.iterations,
            f_reltol=numerics.objective_relative_tolerance,
            g_abstol=numerics.gradient_tolerance / scale,
            show_trace=false,
        ),
    )
    z = project(Optim.minimizer(result))
    gradient = ForwardDiff.gradient(objective, z)
    projected_gradient = norm(z .- clamp.(z .- gradient,
        space.coordinate_lower, space.coordinate_upper), Inf)
    return (
        initial_coordinates=initial, coordinates=z, chi2=objective(z), converged=Optim.converged(result),
        iterations=Optim.iterations(result), projected_gradient=projected_gradient,
    )
end

# ╔═╡ 986dd5d2-bcf8-5301-9191-0313ecc1bf6d
function fit_selected_background(data, space;
                        initials=background_initial_parameters,
                        numerics=default_numerics())
    counts = map(channel -> length(channel.flux), data)
    npoints = sum(counts)
    nfree = length(space.free_indices)
    npoints > nfree ||
        throw(ArgumentError("The selected window needs more points than free parameters"))
    isempty(initials) && throw(ArgumentError("Provide at least one starting point"))
    runs = [fit_background_from(initial, data, space, numerics) for initial ∈ initials]
    finite_indices = findall(run -> isfinite(run.chi2) && all(isfinite, run.coordinates), runs)
    isempty(finite_indices) && error("No finite background fit was obtained")
    converged_indices = filter(i -> runs[i].converged, finite_indices)
    candidates = isempty(converged_indices) ? finite_indices : converged_indices
    best_index = candidates[argmin([runs[i].chi2 for i ∈ candidates])]
    best = runs[best_index]
    predictions = background_predictions(best.coordinates, data, space)
    pulls = background_pulls(best.coordinates, data, space)
    chi2 = map(pull -> sum(abs2, pull), pulls)
    return (
        data=data, counts=counts, npoints=npoints, nfree=nfree, dof=npoints - nfree,
        coordinates=best.coordinates,
        parameter_values=background_parameter_values(best.coordinates, space),
        physical_parameters=physical_background_parameters(best.coordinates, space),
        parameters=decode_background_parameters(best.coordinates, space),
        predictions=predictions, pulls=pulls, chi2=chi2,
        total_chi2=sum(chi2), reduced_chi2=sum(chi2) / (npoints - nfree),
        converged=best.converged, runs=runs, best_run_index=best_index,
    )
end

# ╔═╡ 41e71153-02cd-520b-beb2-fc21ce3a811b
md"""
## 6. Parameter covariance and spectrum uncertainty

The residual Jacobian gives the local joint covariance, converted from optimizer
coordinates to the parameter values in the units listed in the parameter table.
Entry ``(i,j)`` carries the product of those two units. Fixed parameters have zero
covariance rows and columns. The experimental errors set the covariance scale.

`covariance_status` records convergence, rank, conditioning and active bounds.
Parameter errors and optional pointwise spectrum bands are shown when that local
approximation is usable. A rank-deficient fit has `covariance=nothing`.
"""

# ╔═╡ 4baff3dc-6d33-5088-8c6c-675f1c7bfd71
function background_covariance(fit, space; numerics=default_numerics())
    z = fit.coordinates
    nfree, nparameters = length(z), length(space.names)
    if nfree == 0
        return (matrix=zeros(nparameters, nparameters), status=:ok)
    end
    jacobian = ForwardDiff.jacobian(
        coordinates -> background_standardized_residuals(coordinates, fit.data, space), z)
    decomposition = svd(jacobian)
    singular_values = decomposition.S
    threshold = max(numerics.svd_relative_tolerance * maximum(singular_values), eps(Float64))
    count(>(threshold), singular_values) == nfree ||
        return (matrix=nothing, status=:rank_deficient)

    covariance_z = decomposition.V * Diagonal(inv.(singular_values .^ 2)) * decomposition.V'
    transform = ForwardDiff.jacobian(
        coordinates -> background_parameter_values(coordinates, space), z)
    covariance = Matrix(Symmetric(transform * covariance_z * transform'))
    distances = min.(z .- space.coordinate_lower, space.coordinate_upper .- z)
    at_bound = any(distances .<= numerics.bound_tolerance .* (1 .+ abs.(z)))
    condition_number = maximum(singular_values) / minimum(singular_values)
    status = !fit.converged ? :not_converged : at_bound ? :at_bound :
        condition_number >= numerics.max_condition_number ? :ill_conditioned : :ok
    return (matrix=covariance, status=status)
end

# ╔═╡ 931453a6-2122-5afc-8b16-04c39b93ce84
function background_spectrum_band(physical_parameters, covariance, energies, channel)
    parameters = unpack_background_parameters(physical_parameters)
    central = [background_flux(E, parameters, channel) for E ∈ energies]
    if isnothing(covariance)
        return (central=central, standard_error=nothing, lower=nothing, upper=nothing)
    end
    parameter_values = background_parameter_values(physical_parameters)
    standard_error = map(energies) do E
        gradient = ForwardDiff.gradient(parameter_values) do values
            parameters = unpack_background_parameters(physical_background_parameters(values))
            flux_value(background_flux(E, parameters, channel))
        end
        sqrt(max(dot(gradient, covariance * gradient), 0.0)) * differential_flux_unit
    end
    return (central=central, standard_error=standard_error,
            lower=central .- standard_error, upper=central .+ standard_error)
end

# ╔═╡ c7f20d58-be0b-5916-9491-2bee3a98b59d
function background_spectra(physical_parameters, covariance, energies)
    spectra = NamedTuple{background_channels}(Tuple(
        background_spectrum_band(physical_parameters, covariance, energies, channel)
        for channel ∈ background_channels
    ))
    return merge((energy=collect(energies),), spectra)
end

# ╔═╡ b3b627e2-034b-59e1-92a4-8b5e01c85dab
function background_parameter_rows(fit, covariance, space)
    return [
        (
            parameter=name, unit_str=space.metadata[name].unit_str,
            best_fit=fit.parameter_values[i],
            standard_error=haskey(space.fixed, name) ? 0.0 :
                covariance.status == :ok ? sqrt(max(covariance.matrix[i, i], 0.0)) : missing,
            fixed=haskey(space.fixed, name),
        ) for (i, name) ∈ enumerate(space.names)
    ]
end

# ╔═╡ 2d16ac6f-06fc-5e1a-80ca-e34ddc649a09
md"""
## 7. Parameter table, fitted spectra and residuals

`plot_background_fit(result)` draws an existing fit without rerunning the optimizer.
Parameters and prepared data can also be supplied directly:
```julia
plot_background_fit(result.parameters, result.data;
    E_min=result.E_min,
    covariance=result.covariance)
```
For this direct form, the optional covariance uses the parameter order and units
in `background_parameter_metadata`. Without covariance, the central curves and
data residuals are still drawn. `bands=false` hides the model uncertainty bands.

The spectrum panel compares the model with all retained bins above `E_min`,
with each curve limited to the measured energy range of its data channel.
Its shading is the pointwise local ``1\sigma`` model uncertainty.
The residual panel shows
``r_i=(J_i-F_i)/\sigma_i`` with data error bars of size one and a grey
``\pm1`` reference strip. Coloured bands around zero show
``\pm\sigma_F(E_i)/\sigma_i`` at the data energies.
These model bands are displayed separately from the data errors and are not
confidence intervals for the fitted residuals.
"""

# ╔═╡ 39b7f345-0f20-5692-944c-d491b90c175a
function parameter_table(rows)
    number(value) = ismissing(value) ? "unavailable" : string(round(value; sigdigits=5))
    header = "| Parameter | Best fit | Local 1σ | Unit |\n|---|---:|---:|---|\n"
    lines = ["| `$(row.parameter)` | $(number(row.best_fit)) | " *
             "$(row.fixed ? "fixed" : number(row.standard_error)) | $(row.unit_str) |"
             for row ∈ rows]
    return Markdown.parse(header * join(lines, '\n'))
end

# ╔═╡ 0af5046c-7b1a-428a-9d58-8d8d41533f1a
function background_residuals(physical_parameters, data, covariance=nothing)
    function dimensionless(value)
        @check_EU_dimension value 0
        return EUval(value)
    end
    return NamedTuple{keys(data)}(map(keys(data)) do channel
        measurements = data[channel]
        prediction = background_spectrum_band(
            physical_parameters, covariance, measurements.energy, channel)
        central = dimensionless.((measurements.flux .- prediction.central) ./ measurements.uncertainty)
        model_error = isnothing(prediction.standard_error) ? nothing :
            dimensionless.(prediction.standard_error ./ measurements.uncertainty)
        (central=central, model_standard_error=model_error)
    end)
end

# ╔═╡ c85500cc-54f9-5527-97ec-e211e6565a75
function make_background_fit_figure(data, residuals, spectra;
                                    summary="Supplied background parameters")
    figure = Figure(size=(1050, 760), fontsize=17)
    channel_labels = (electron="AMS e−", positron="AMS e+", electron_positron="DAMPE e− + e+")
    channel_summary = join([
        "$(channel_labels[channel]): χ²=$(round(sum(abs2, residuals[channel].central); digits=2)), N=$(length(data[channel].flux))"
        for channel ∈ keys(data)], "   |   ")
    has_bands = !isnothing(spectra.electron.standard_error)
    Label(figure[0, 1], summary * "\n" * channel_summary * "\n" *
          (has_bands ? "Spectrum shading: local ±1σ model uncertainty." : "Model uncertainty bands not shown.");
          fontsize=14, tellwidth=false)
    flux_label = L"E^3 J\;[\mathrm{GeV}^2\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}\,\mathrm{sr}^{-1}]"
    flux_axis = Axis(figure[1, 1]; title="Fit and retained measurements",
                    ylabel=flux_label, xscale=log10, yscale=log10)
    pull_title = has_bands ? "Grey: ±1 reference. Coloured: model ±σ / data σ." : "Grey: ±1 reference."
    pull_axis = Axis(figure[2, 1]; title=pull_title, titlesize=13,
                     xlabel=L"E\;[\mathrm{GeV}]",
                     ylabel=L"(J_\mathrm{data}-J_\mathrm{model})/\sigma", xscale=log10)
    colors = (electron=:seagreen4, positron=:darkorange2, electron_positron=:mediumpurple4)
    labels = (electron="e−", positron="e+", electron_positron="e− + e+")
    markers = (electron=:circle, positron=:utriangle, electron_positron=:diamond)
    energies = EUval.(Ref(GeV), spectra.energy)
    for channel ∈ keys(data)
        prediction, color = spectra[channel], colors[channel]
        displayed = energies .<= EUval(GeV, maximum(data[channel].energy_max))
        channel_energies = energies[displayed]
        if !isnothing(prediction.standard_error)
            lower = flux_value.(prediction.lower[displayed])
            lower = ifelse.(lower .> 0, lower, NaN)
            band!(flux_axis, channel_energies, channel_energies .^ 3 .* lower,
                  channel_energies .^ 3 .* flux_value.(prediction.upper[displayed]); color=(color, 0.17))
        end
        lines!(flux_axis, channel_energies, channel_energies .^ 3 .* flux_value.(prediction.central[displayed]);
               color=color, linewidth=2, label=labels[channel])
    end
    hspan!(pull_axis, [-1.0], [1.0]; color=(:grey, 0.12))
    for channel ∈ keys(data)
        measurements, color = data[channel], colors[channel]
        data_energies = EUval.(Ref(GeV), measurements.energy)
        residual = residuals[channel]
        if !isnothing(residual.model_standard_error)
            band!(pull_axis, data_energies, -residual.model_standard_error,
                  residual.model_standard_error; color=(color, 0.17))
        end
        errorbars!(flux_axis, data_energies, data_energies .^ 3 .* flux_value.(measurements.flux),
                   data_energies .^ 3 .* flux_value.(measurements.uncertainty); color=(color, 0.65))
        scatter!(flux_axis, data_energies, data_energies .^ 3 .* flux_value.(measurements.flux);
                 color=color, marker=markers[channel], markersize=6)
        errorbars!(pull_axis, data_energies, residual.central, ones(length(data_energies));
                   color=(color, 0.55), whiskerwidth=3)
        scatter!(pull_axis, data_energies, residual.central;
                 color=color, marker=markers[channel], markersize=6)
    end
    hlines!(pull_axis, [0.0]; color=(:black, 0.6), linestyle=:dash)
    for axis ∈ (flux_axis, pull_axis)
        xlims!(axis, first(energies), last(energies))
    end
    linkxaxes!(flux_axis, pull_axis)
    hidexdecorations!(flux_axis; grid=false)
    axislegend(flux_axis; position=:lb)
    rowsize!(figure.layout, 1, Relative(0.65))
    return figure
end

# ╔═╡ 490fbbe5-ec8d-4d97-a037-09e85221f88a
begin
    """
        plot_background_fit(parameters, data; E_min, covariance=nothing, ...)
        plot_background_fit(result; bands=true, spectrum_points=160)

    Compare supplied physical parameters with prepared measurements, or plot a
    `fit_background` result. No optimization is performed. Return a Makie Figure
    with fitted spectra and standardized residuals. An optional numeric
    covariance in the order and units of `background_parameter_metadata` supplies
    local model uncertainty bands. The result form uses its covariance only when
    `covariance_status == :ok`.
    """
    function plot_background_fit(parameters::NamedTuple, data::NamedTuple;
                                 E_min, covariance=nothing, bands=true,
                                 spectrum_points=160,
                                 summary="Supplied background parameters")
        validate_background_threshold(E_min)
        selected = select_background_data(data, E_min)
        band_covariance = bands ? covariance : nothing
        spectra = background_spectra(parameters, band_covariance,
                                    background_energy_grid(selected, E_min; points=spectrum_points))
        residuals = background_residuals(parameters, selected, band_covariance)
        return make_background_fit_figure(selected, residuals, spectra; summary)
    end

    function plot_background_fit(result; bands=true, spectrum_points=160)
        covariance = result.covariance_status == :ok ? result.covariance : nothing
        summary = "$(result.strategy)   E_min = $(EUval(GeV, result.E_min)) GeV. " *
                  "χ²/dof = $(round(result.reduced_chi2; sigdigits=5)). " *
                  "Converged: $(result.converged). Covariance: $(result.covariance_status)."
        return plot_background_fit(result.parameters, result.data;
            E_min=result.E_min,
            covariance, bands, spectrum_points, summary)
    end
end

# ╔═╡ 302c68ac-2fd9-54e2-92d7-1f6fc4fb67ad
md"""
## 8. One-call interface

The lower threshold `E_min` is required. Each call fits its own parameters and
returns `reduced_chi2`, `parameters`, `covariance`, `parameter_table`, and `figure`.
The covariance row/column order and units are returned as `parameter_names`
and `parameter_units`. `bands=false` hides spectrum error bands.

Example calls with a 50 GeV lower threshold:
```julia
AMS_fit = fit_background(strategy=:AMS, E_min=GeV(50));
DAMPE_fit = fit_background(strategy=:AMS_DAMPE, E_min=GeV(50));
plot_background_fit(AMS_fit)
plot_background_fit(DAMPE_fit)
```
"""

# ╔═╡ d6effab7-80e6-57ff-aaee-028806ecdd65
"""
    fit_background(; E_min, strategy=:AMS, ...)

Fit all complete bins above E_min using AMS electrons + positrons (`:AMS`) or
AMS positrons + DAMPE total flux (`:AMS_DAMPE`). Both charge parameter sets
are fitted together. Optional `bounds`, `fixed`, `initials`, and `numerics`
override the settings above. `data` may supply prepared channels for the fit.

Return best-fit physical parameters, chi-square summaries, local covariance,
a parameter table, spectra and a figure with standardized residuals.
Covariance entries use the products of the returned parameter units.
`make_figure=false` skips figure construction when only numerical output is needed.
"""
function fit_background(; E_min, strategy=:AMS, data=nothing,
                        bounds=background_parameter_bounds,
                        fixed=background_fixed_parameters,
                        initials=background_initial_parameters,
                        numerics=default_numerics(), bands=true,
                        spectrum_points=160, make_figure=true)
    strategy ∈ keys(background_fit_strategies) ||
        throw(ArgumentError("Use strategy=:AMS or :AMS_DAMPE"))
    validate_background_threshold(E_min)
    channels = background_fit_strategies[strategy]
    inputs = isnothing(data) ? load_background_data(strategy) : NamedTuple{channels}(data)
    selected = select_background_data(inputs, E_min)
    for channel ∈ channels
        isempty(selected[channel].flux) &&
            throw(ArgumentError("No complete bins for $channel in the specified window"))
    end
    space = background_parameter_space(bounds, fixed)
    fit = fit_selected_background(selected, space; initials, numerics)
    covariance = background_covariance(fit, space; numerics)
    band_covariance = bands && covariance.status == :ok ? covariance.matrix : nothing
    spectra = background_spectra(fit.physical_parameters, band_covariance,
                                background_energy_grid(selected, E_min; points=spectrum_points))
    rows = background_parameter_rows(fit, covariance, space)
    result = (
        strategy=strategy, E_min=E_min, parameters=fit.physical_parameters,
        chi2=fit.total_chi2, reduced_chi2=fit.reduced_chi2, dof=fit.dof,
        npoints=fit.npoints, nfree=fit.nfree, counts=fit.counts, channel_chi2=fit.chi2,
        converged=fit.converged, covariance=covariance.matrix,
        covariance_status=covariance.status, parameter_names=space.names,
        parameter_units=map(entry -> entry.unit_str, space.metadata),
        parameter_table=parameter_table(rows), parameter_rows=rows,
        data=selected, pulls=fit.pulls, spectra=spectra,
    )
    figure = make_figure ? plot_background_fit(result; bands, spectrum_points) : nothing
    return merge(result, (figure=figure,))
end

# ╔═╡ 8e04c129-6164-563c-9409-83d01056c025
# Run configured analyses.
md"""
## 9. Run a specified lower threshold

Set `E_min` below to a NaturalUnits energy and choose a strategy.
The figures and covariance refer to all selected data above that threshold. The full-data visual
comparison is in `data-plot.pluto.jl`.
"""

# ╔═╡ 603152f8-64ab-5351-8e32-e638194de6b4
background_fit_settings = (
    strategy=:AMS_DAMPE,
    E_min=GeV(40), # Use GeV(40) for a lower-threshold comparison.
)

# ╔═╡ d06b0ff8-e6fb-551b-af65-7a1b237ebdd1
background_fit_result = if isnothing(background_fit_settings.E_min)
    nothing
else
    fit_background(; background_fit_settings...)
end;

# ╔═╡ c829080a-365b-563c-8e0a-9023f2dfbe8a
if isnothing(background_fit_result)
    md"Set `E_min` above to run the fit."
else
    result = background_fit_result
    Markdown.parse("""
    **$(result.strategy)**: χ² = $(result.chi2), dof = $(result.dof),
    χ²/dof = $(result.reduced_chi2). Converged: $(result.converged).
    Covariance: `$(result.covariance_status)`.
    """)
end

# ╔═╡ 6a7f2c18-ab5d-5567-813a-3e766ade585f
isnothing(background_fit_result) ? nothing : background_fit_result.parameter_table

# ╔═╡ f6e77c35-ca37-5286-9e1b-c7340324c80f
isnothing(background_fit_result) ? nothing : background_fit_result.figure

# ╔═╡ c51fdffa-c1ab-5007-8514-5c5fa96f5bb3
md"Covariance rows and columns follow the parameter table. Each entry has the product of the two listed units."

# ╔═╡ 5bdbb224-505a-5379-a1df-c70dccef3ebf
isnothing(background_fit_result) ? nothing : background_fit_result.covariance

# ╔═╡ Cell order:
# ╠═97dd5ef4-6bb5-5531-829d-7352f2044073
# ╠═c2162faf-355b-57a8-9d79-52e62f03d31c
# ╟─1955b0f1-b919-518c-8c5b-43e2f7b5a6dc
# ╟─0ed3e31c-7a52-56d2-98b8-ef905086afe3
# ╠═41d8f99a-f021-5112-9158-4036bc61660a
# ╠═3c6ab0db-9a1c-558d-8338-c1eeb71b1d0b
# ╠═b5ce58d3-1cd3-53dd-9a48-937949dea816
# ╟─92fdfccd-5471-591f-b2db-60788043e605
# ╟─3744a5cc-61e6-4cc5-a2a9-1bac68a9cb63
# ╠═8f0235b8-7f17-5e00-9941-e93fb48933b2
# ╠═58f03616-ac8d-5c17-822c-5556a8f32ca7
# ╠═ff397277-672e-5813-8fe0-a4d6ebc1606e
# ╠═4ff905dc-f2ef-5ec2-ab25-59bb5f483e51
# ╠═933450cd-0ee2-5fc2-add2-63ebe9da9efa
# ╟─3dd95f9b-4f2c-56cf-9adf-735bee6ea300
# ╠═598bba1a-953b-5085-93aa-042c4366f6fb
# ╠═01cd9677-d31a-5cc9-849f-f89983ff1668
# ╠═67e02064-abbe-51c7-bbc9-af8ee56b45d7
# ╠═d15f2c11-32bb-59c0-908c-35d9ee1e3a44
# ╠═ac32361e-2866-5e20-8d43-907b8a049b2e
# ╠═fcd83373-6b69-576b-8eff-f440fd700da9
# ╠═9d5b5516-2a36-5f28-a6a6-12d4dee1e6f3
# ╟─fef949b0-4cad-54b6-8dd0-4354032c2919
# ╠═1f3cbf9f-7420-5569-a632-2aaf97ec9b19
# ╠═c2defc4c-c686-57a7-acba-4c217c367bc5
# ╠═aeb77dc9-f168-5529-b9ee-eaaced92bf64
# ╠═484285e8-90e0-526b-94ef-47d81e21c674
# ╠═ce6d9c40-aedd-5ad5-8595-818b8db807c5
# ╠═a936d1ed-3d0e-5a53-a934-1c99fe6a469b
# ╠═d36baa99-5336-5c88-b224-471bfd5ea85d
# ╠═2fc9db5e-1b73-581e-8905-5a1a09be12f4
# ╠═57e0e552-d709-5a89-9939-ebddd38a46c0
# ╠═c271919b-f8a3-5961-aef4-8fc8f854c0d8
# ╟─e9deb9f2-00c3-504c-834f-7655fa61be56
# ╠═aeae588a-599b-549a-9deb-39b7fd147d4b
# ╠═986dd5d2-bcf8-5301-9191-0313ecc1bf6d
# ╟─41e71153-02cd-520b-beb2-fc21ce3a811b
# ╠═4baff3dc-6d33-5088-8c6c-675f1c7bfd71
# ╠═931453a6-2122-5afc-8b16-04c39b93ce84
# ╠═c7f20d58-be0b-5916-9491-2bee3a98b59d
# ╠═b3b627e2-034b-59e1-92a4-8b5e01c85dab
# ╟─2d16ac6f-06fc-5e1a-80ca-e34ddc649a09
# ╠═39b7f345-0f20-5692-944c-d491b90c175a
# ╠═0af5046c-7b1a-428a-9d58-8d8d41533f1a
# ╠═c85500cc-54f9-5527-97ec-e211e6565a75
# ╠═490fbbe5-ec8d-4d97-a037-09e85221f88a
# ╟─302c68ac-2fd9-54e2-92d7-1f6fc4fb67ad
# ╠═d6effab7-80e6-57ff-aaee-028806ecdd65
# ╟─8e04c129-6164-563c-9409-83d01056c025
# ╟─603152f8-64ab-5351-8e32-e638194de6b4
# ╟─d06b0ff8-e6fb-551b-af65-7a1b237ebdd1
# ╟─c829080a-365b-563c-8e0a-9023f2dfbe8a
# ╟─6a7f2c18-ab5d-5567-813a-3e766ade585f
# ╠═f6e77c35-ca37-5286-9e1b-c7340324c80f
# ╟─c51fdffa-c1ab-5007-8514-5c5fa96f5bb3
# ╠═5bdbb224-505a-5379-a1df-c70dccef3ebf
