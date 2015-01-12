# Generate the modules.
include("JFFoundationModule.jl")
include("PhysicalUnitModule.jl")
include("RotationUtilModule.jl")
include("FESetModule.jl")
include("FENodeSetModule.jl")
include("NodalFieldModule.jl")
include("MeshUtilModule.jl")
include("MeshSelectionModule.jl")
include("MeshExportModule.jl")
include("MeshModificationModule.jl")
include("MeshTransformModule.jl")
include("MeshQuadrilateralModule.jl")
include("MeshLineModule.jl")
include("MeshTriangleModule.jl")
include("MeshHexahedronModule.jl")
include("ForceIntensityModule.jl")
include("PropertyHeatDiffusionModule.jl")
include("MaterialHeatDiffusionModule.jl")
include("PropertyAcousticFluidModule.jl")
include("MaterialAcousticFluidModule.jl")
include("AssemblyModule.jl")
include("IntegRuleModule.jl")
include("FEMMBaseModule.jl")
include("FEMMHeatDiffusionModule.jl")
include("FEMMAcousticsModule.jl")
include("HeatDiffusionAlgorithmModule.jl")
include("AcousticsAlgorithmModule.jl")

"""

JFinEALE (C) 2014-2015, Petr Krysl
  
"""

# This module is defined simply in order to make the package manager happy.
module JFinEALE
end


