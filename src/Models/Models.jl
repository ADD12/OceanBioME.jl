module Models

export Sediments

export NPZD, 
       NutrientPhytoplanktonZooplanktonDetritus, 
       LOBSTER,
       PISCES

export SLatissima

export CarbonChemistry

export GasExchange, 
       CarbonDioxideGasExchangeBoundaryCondition, 
       OxygenGasExchangeBoundaryCondition, 
       GasExchangeBoundaryCondition,
       ScaledGasTransferVelocity,
       SchmidtScaledTransferVelocity,
       CarbonDioxidePolynomialSchmidtNumber,
       OxygenPolynomialSchmidtNumber

include("Sediments/Sediments.jl")
include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SLatissima.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")

using .Sediments
using .LOBSTERModel
using .NPZDModel
using .PISCESModel
using .SLatissimaModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
