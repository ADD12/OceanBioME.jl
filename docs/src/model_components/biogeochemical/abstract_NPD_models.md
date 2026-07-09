# [Nutrients-Plankton-Detritus framework](@id npd_framework)

The Nutrient-Plankton-Detritus framework allows biogeochemical models to be built from shared components. Individual biogeochemical models like [LOBSTER](@ref LOBSTER) or [NPZD](@ref NPZD) can be thought of as 'recipes' that combine these ingredients in different ways. 

The framework includes five component types: Nutrients, Plankton (including phytoplankton and zooplankton), Detritus, Inorganic Carbon, and Oxygen. This framework is extremely flexible. Models built using this framework don't need to use all of these component types, but they can also include mutiple components. For example, multiple nutrient and phytoplankton types can be used in the biogeochemical model and it is even possible not to include explicit phytoplankton (as in [ImplicitBiology](@ref ImplicitBiology)). This framework allows biogeochemical models to be highly composable, and greatly aids development of new biogeochemical models. For a step-by-step guide to creating a new biogeochemical model using the Nutrient-Plankton-Detritus framework, see [Implementing new models](@ref model_implementation).

The `NutrientsPlanktonDetritus` framework organises the biogeochemistry into five component types. Biogeochemical models must specify which of the existing types are used for each component (or use the default of `nothing` in the case of Inorganic Carbon and Oxygen).

| Component | Keyword | Available types |
|------|---------|-----------------|
| Nutrients | `nutrients` | `Nutrients` (constructed with `NitrateAmmonia`, `N`, `PO₄`, `Fe`, or `Si` sub-components) |
| Plankton | `plankton` | `Abiotic`, `ImplicitProductivity`, `PhytoZoo` |
| Detritus | `detritus` | `InstantRemineralisationDetritus`, `Detritus`, `DissolvedParticulate`, `CarbonNitrogenDissolvedParticulate` |
| Inorganic Carbon | `inorganic_carbon` | `nothing` (default), `CarbonateSystem` |
| Oxygen | `oxygen` | `nothing` (default), `Oxygen` |

## Preset constructors

For common configurations, convenience constructors assemble the components using the 'recipe' for that biogeochemical model with default parameters. You can build each biogeochemical model with default parameters simply by passing the model grid as an argument:

- `LOBSTER(grid)` — medium-complexity model with phytoplankton, zooplankton, nitrate, ammonia, and dissolved/particulate detritus. See the [LOBSTER](@ref LOBSTER) page.
- `NPZD(grid)` — simple four-compartment nutrient–phytoplankton–zooplankton–detritus model from [Kuhn2015](@citet). See the [NPZD](@ref NPZD) page.
- `ImplicitBiology(grid)` — nutrient-limited community productivity with no explicit plankton biomass. See the [ImplicitBiology](@ref ImplicitBiology) page.

Users can over-ride the default parameters by optionally passing in parameter values. Not all parameters need to be specified.  Any parameters that are not specified will be set to the default value. For example:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid;
                            biogeochemistry = NPZD(grid;
                                                   plankton = PhytoZoo(grid;
                                                                      phytoplankton_maximum_growth_rate = 1 / day,
                                                                      maximum_grazing_rate = 2 / day)))
```

All preset constructors accept the same optional keyword arguments as `NutrientsPlanktonDetritus` itself — in particular `inorganic_carbon`, `oxygen`, `light_attenuation`, `sediment`, and `scale_negatives` — so you can extend any of them without building from scratch. For example, the following code adds inorganic carbon and oxygen to the LOBSTER model

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid;
                            biogeochemistry = LOBSTER(grid;
                                                      inorganic_carbon = CarbonateSystem(),
                                                      oxygen = Oxygen()))
```

## Custom configurations

You can create a fully custom biogeochemical model by calling `NutrientsPlanktonDetritus` directly and specifying each component. For example, to build a model with Nitrate, Ammonia, Phosphate, and Iron, one phytoplankton and one zooplankton type, one dissolved detritus pool, and an active carbonate system: 

```julia
biogeochemistry = NutrientsPlanktonDetritus(grid;
                                            nutrients  = Nutrients(NitrateAmmonia(), PO₄, Fe, nothing),
                                            plankton   = PhytoZoo(grid),
                                            detritus   = DissolvedParticulate(grid),
                                            inorganic_carbon = CarbonateSystem())
```

The framework dispatches through the components to assemble the required tracers and auxiliary fields automatically, so all standard Oceananigans boundary conditions and output writers work as normal.

