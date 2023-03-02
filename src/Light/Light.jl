"
Light attenuation by chlorophyll as described by [Karleskind2011](@cite) (implimented as twoBand) and [Morel1988](@cite).

These should be used by setting up a callback like:
`simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)), TimeStepCallsite())`
"
module Light
export TwoBandPhotosyntheticallyActiveRatiation, update_PAR!

using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans: Center, Face, fields
using Oceananigans.Grids: xnode, ynode, znodes
using Oceananigans.Fields: CenterField, TracerFields
using Oceananigans.BoundaryConditions: fill_halo_regions!, 
                                       ValueBoundaryCondition, 
                                       FieldBoundaryConditions, 
                                       regularize_field_boundary_conditions, 
                                       ContinuousBoundaryFunction
using Oceananigans.Units

import Adapt: adapt_structure
import Base: show, summary

import Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields
import Oceananigans.BoundaryConditions: _fill_top_halo!

# Fallback
update_PAR!(model, PAR) = nothing
required_PAR_fields(PAR) = ()

include("2band.jl")
include("morel.jl")

end
