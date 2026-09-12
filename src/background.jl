# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

################################################################################
export background_flux

"""
    background_flux(energy::EnergyUnit, param::NamedTuple; reference=GeV(1.0))

Evaluate the empirical broken-power-law differential flux

```math
J(E) = C (E/E_0)^{-\\gamma}
        \\left[1 + (E/E_\\mathrm{br})^{\\Delta\\gamma}\\right]^{-1}.
```

`param` contains `C`, `Ebr`, `γ`, and `Δγ`. `C` has natural-unit mass
dimension +2. `energy`, `Ebr`, and `reference` have dimension +1, and the
indices are dimensionless real numbers. All values must be finite, with
positive `C`, `Δγ`, and `reference`, and total energies `energy` and `Ebr`
above the electron rest energy.

`reference` is the normalization energy `E₀`. Solar modulation is not included.
The returned flux uses the canonical [`EU`](@ref) basis. Fit thresholds and
numerical search bounds are specified by the caller.
"""
function background_flux(energy::EnergyUnit, param::NamedTuple;
                            reference::EnergyUnit=GeV(1.0))
    E = _canonical_unit(energy, 1, "total energy")
    Ebr = _canonical_unit(param.Ebr, 1, "break energy")
    E0 = _canonical_unit(reference, 1, "reference energy")
    C = _canonical_unit(param.C, 2, "background normalization")
    _require_dimension(param.γ, 0, "spectral index")
    _require_dimension(param.Δγ, 0, "spectral index change")
    param.γ isa Real && param.Δγ isa Real ||
        throw(ArgumentError("Spectral indices must be real numbers"))

    _require_finite(E, "total energy")
    _require_finite(Ebr, "break energy")
    rest_energy = convert(EU, electron_mass)
    E > rest_energy ||
        throw(DomainError(energy, "Total energy must exceed the rest energy"))
    Ebr > rest_energy || throw(DomainError(param.Ebr,
                                "Break energy must exceed the rest energy"))
    _require_positive_finite(E0, "reference energy")
    _require_positive_finite(C, "background normalization")
    _require_finite(param.γ, "spectral index")
    _require_positive_finite(param.Δγ, "spectral index change")

    # Stable log(1 + exp(x)), including AD through the spectral break.
    x = param.Δγ * log(EUval(E / Ebr))
    softplus = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
    return C * exp(-param.γ * log(EUval(E / E0)) - softplus)
end
################################################################################
