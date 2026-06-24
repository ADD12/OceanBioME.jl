# [Implementing new models](@id model_implementation)

There are two main ways to extend OceanBioME with new biology:

1. **Add a new plankton component to the NPD framework** — implement a `plankton` type that slots into the existing [Nutrients-Plankton-Detritus framework](@ref npd_framework). The framework automatically handles nutrient uptake, detritus production, inorganic carbon, and oxygen coupling.

2. **Implement a completely new BGC model** — subtype `AbstractContinuousFormBiogeochemistry` from Oceananigans for full control over all tracer tendencies. This is more work but imposes no constraints on model structure.

This page focuses on the first approach, which is appropriate for most new plankton models.

## The NPD plankton interface

A plankton component must implement the following four functions, all called with signature `(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)`:

| Function | Returns | Description |
|----------|---------|-------------|
| `nutrient_uptake` | mmol N / m³ / s | Total N removed from the inorganic nutrient pool by growth |
| `dissolved_waste` | mmol N / m³ / s | N exported to dissolved organic pools (exudate, dissolved mortality) |
| `solid_waste` | mmol N / m³ / s | N exported to particulate organic pools (mortality, fecal pellets) |
| `inorganic_waste` | mmol N / m³ / s | N returned directly to the inorganic pool (excretion, respiration) |

The `nutrient_uptake` function may optionally be specialised per tracer by adding a `::Val{:NO₃}` (or `::Val{:NH₄}` etc.) argument between `grid` and `plankton`, which is useful when ammonia and nitrate are tracked separately.

Additionally:

- `required_biogeochemical_tracers` must return a tuple of the tracer names the component owns (e.g. `(:P, :Z)`).
- `required_biogeochemical_auxiliary_fields` must return a tuple of auxiliary fields needed (typically `(:PAR,)`).
- A tracer tendency method must be defined for each owned tracer via `(bgc::NutrientsPlanktonDetritus)(i, j, k, grid, ::Val{:P}, ...)`.

The default elemental ratios (Redfield: C:N:P:Fe = 106:16:1:0.0032) are used automatically unless you override `carbon_ratio`, `nitrogen_ratio`, `phosphate_ratio`, or `iron_ratio`. You may also define `detritus_grazing` to implement zooplankton-like grazing on the detritus pools.

## Example: simple phytoplankton

Here we implement a minimal phytoplankton model with Michaelis-Menten light and nutrient limitation, linear mortality, and exudation of dissolved organic matter.

### Imports

```julia
using OceanBioME, Oceananigans
using Oceananigans.Units
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields

# Import the NPD interface functions we are adding methods to
import OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels: inorganic_waste, nutrient_uptake
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: dissolved_waste, solid_waste
```

### The struct

```julia
@kwdef struct SimplePhytoplankton{FT}
    maximum_growth_rate      :: FT = 1.5 / day   # 1/s
    light_half_saturation    :: FT = 30.0         # W/m²
    nitrate_half_saturation  :: FT = 0.5          # mmol N/m³
    mortality_rate           :: FT = 0.01 / day   # 1/s
    exudate_fraction         :: FT = 0.1          # fraction of gross growth exuded as DOM
end

required_biogeochemical_tracers(::SimplePhytoplankton)        = (:P,)
required_biogeochemical_auxiliary_fields(::SimplePhytoplankton) = (:PAR,)
```

### Growth rate

We compute gross phytoplankton growth separately so it can be reused:

```julia
@inline function gross_growth(i, j, k, grid, p::SimplePhytoplankton, fields, auxiliary_fields)
    PAR = @inbounds auxiliary_fields.PAR[i, j, k]
    N   = @inbounds fields.NO₃[i, j, k]
    P   = @inbounds fields.P[i, j, k]

    L_light    = PAR / (PAR + p.light_half_saturation)
    L_nutrient = N   / (N   + p.nitrate_half_saturation)

    return p.maximum_growth_rate * L_light * L_nutrient * P
end
```

### Tracer tendency

```julia
@inline function (bgc::NutrientsPlanktonDetritus)(i, j, k, grid, ::Val{:P}, clock, fields, auxiliary_fields)
    p        = bgc.plankton
    growth   = gross_growth(i, j, k, grid, p, fields, auxiliary_fields)
    exudate  = p.exudate_fraction * growth
    mortality = p.mortality_rate * @inbounds fields.P[i, j, k]

    return growth - exudate - mortality
end
```

### Interface methods

```julia
# Nutrient uptake: all gross growth removes N from the nutrient pool
@inline nutrient_uptake(i, j, k, grid, p::SimplePhytoplankton, bgc, fields, aux) =
    gross_growth(i, j, k, grid, p, fields, aux)

# Specialise by tracer when NitrateAmmonia nutrients are used (no ammonia uptake here)
@inline nutrient_uptake(i, j, k, grid, ::Val{:NO₃}, p::SimplePhytoplankton, bgc, fields, aux) =
    gross_growth(i, j, k, grid, p, fields, aux)
@inline nutrient_uptake(i, j, k, grid, ::Val{:NH₄}, p::SimplePhytoplankton, bgc::NutrientsPlanktonDetritus{FT}, fields, aux) where FT =
    zero(FT)

# Dissolved waste: exudate from growth
@inline dissolved_waste(i, j, k, grid, p::SimplePhytoplankton, bgc, fields, aux) =
    p.exudate_fraction * gross_growth(i, j, k, grid, p, fields, aux)

# Solid waste: mortality goes to particulate pool
@inline solid_waste(i, j, k, grid, p::SimplePhytoplankton, bgc, fields, aux) =
    p.mortality_rate * @inbounds fields.P[i, j, k]

# No direct inorganic release
@inline inorganic_waste(i, j, k, grid, ::SimplePhytoplankton, bgc::NutrientsPlanktonDetritus{FT}, args...) where FT =
    zero(FT)
```

The N budget closes: `nutrient_uptake = dP/dt + dissolved_waste + solid_waste + inorganic_waste`.

### Using the new model

The component can be plugged into any NPD preset constructor or into `NutrientsPlanktonDetritus` directly:

```julia
grid = RectilinearGrid(size = (1, 1, 32), extent = (1, 1, 200))

# Drop SimplePhytoplankton into LOBSTER's nutrient/detritus setup
biogeochemistry = NutrientsPlanktonDetritus(grid;
                                            nutrients = Nutrients(NitrateAmmonia(), nothing, nothing, nothing),
                                            plankton  = SimplePhytoplankton(),
                                            detritus  = DissolvedParticulate(grid))

model = NonhydrostaticModel(grid; biogeochemistry)
```

Because the NPD framework dispatches through your interface methods to build all the tracer tendencies, no further coupling code is needed.

## GPU support

To run on a GPU you must tell Julia how to transfer your struct to the device. Add this after your struct definition:

```julia
import Adapt: adapt_structure
using Adapt

Adapt.adapt_structure(to, p::SimplePhytoplankton) =
    SimplePhytoplankton(adapt(to, p.maximum_growth_rate),
                        adapt(to, p.light_half_saturation),
                        adapt(to, p.nitrate_half_saturation),
                        adapt(to, p.mortality_rate),
                        adapt(to, p.exudate_fraction))
```

## Implementing a completely new BGC model

If your model does not fit the NPD component structure — for example if it requires non-standard tracer coupling or has no plankton at all — you can implement it as a standalone `AbstractContinuousFormBiogeochemistry`. See the [Oceananigans biogeochemistry documentation](https://clima.github.io/OceananigansDocumentation/stable/) for the full interface, and the OceanBioME `Biogeochemistry` wrapper to add light attenuation, sediments, and particles.
