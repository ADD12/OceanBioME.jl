using OceanBioME, Test, CUDA, Oceananigans, JLD2, Oceananigans.Units, Documenter

using Oceananigans.Fields: ConstantField

architecture = CUDA.functional() ? GPU() : CPU()