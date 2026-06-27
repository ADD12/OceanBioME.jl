# [Nutrients-Plankton-Detritus framework](@id npd_framework)

[LOBSTER](@ref LOBSTER), [NPZD](@ref NPZD), and [ImplicitBiology](@ref ImplicitBiology) are all built on a common modular framework called `NutrientsPlanktonDetritus`. This framework organises the biogeochemistry into five pluggable component slots:

| Slot | Keyword | Available types |
|------|---------|-----------------|
| Nutrients | `nutrients` | `Nutrients` (constructed with `NitrateAmmonia`, `N`, `PO₄`, `Fe`, or `Si` sub-components) |
| Plankton | `plankton` | `Abiotic`, `ImplicitProductivity`, `PhytoZoo` |
| Detritus | `detritus` | `InstantRemineralisation`, `Detritus`, `DissolvedParticulate`, `VariableDissolvedParticulate` |
| Inorganic carbon | `inorganic_carbon` | `nothing` (default), `CarbonateSystem` |
| Oxygen | `oxygen` | `nothing` (default), `Oxygen` |

## Preset constructors

For common configurations, convenience constructors assemble these components with sensible defaults. You can build each by passing a grid:

- `LOBSTER(grid)` — medium-complexity model with phytoplankton, zooplankton, nitrate, ammonia, and dissolved/particulate detritus. See the [LOBSTER](@ref LOBSTER) page.
- `NPZD(grid)` — simple four-compartment nutrient–phytoplankton–zooplankton–detritus model from [Kuhn2015](@citet). See the [NPZD](@ref NPZD) page.
- `ImplicitBiology(grid)` — nutrient-limited community productivity with no explicit plankton biomass. See the [ImplicitBiology](@ref ImplicitBiology) page.

All preset constructors accept the same optional keyword arguments as `NutrientsPlanktonDetritus` itself — in particular `inorganic_carbon`, `oxygen`, `light_attenuation`, `sediment`, and `scale_negatives` — so you can extend any of them without building from scratch:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid;
                            biogeochemistry = LOBSTER(grid;
                                                      inorganic_carbon = CarbonateSystem(),
                                                      oxygen = Oxygen()))
```

## Custom configurations

You can assemble a fully custom model by calling `NutrientsPlanktonDetritus` directly and specifying each component:

```julia
biogeochemistry = NutrientsPlanktonDetritus(grid;
                                            nutrients  = Nutrients(NitrateAmmonia(), PO₄, Fe, nothing),
                                            plankton   = PhytoZoo(grid),
                                            detritus   = DissolvedParticulate(grid),
                                            inorganic_carbon = CarbonateSystem())
```

The framework dispatches through the components to assemble the required tracers and auxiliary fields automatically, so all standard Oceananigans boundary conditions and output writers work as normal.

For a step-by-step guide to implementing your own plankton component, see [Implementing new models](@ref model_implementation).
