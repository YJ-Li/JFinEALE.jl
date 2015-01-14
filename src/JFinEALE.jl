module JFinEALE
"""

JFinEALE (C) 2014-2015, Petr Krysl
  
"""

# Generate the modules.

include("JFFoundationModule.jl")
using JFinEALE.JFFoundationModule
export JFInt, JFFlt, JFCplxFlt, JFFltVec, JFIntVec, JFFltMat, JFIntMat, JFMat, JFVec
export ModelDataDictionary, dadd!, dmake

include("PhysicalUnitModule.jl")
using JFinEALE.PhysicalUnitModule
phun=JFinEALE.PhysicalUnitModule.phun
export phun

include("RotationUtilModule.jl")
using JFinEALE.RotationUtilModule
export rotmat, skewmat!, cross

include("FESetModule.jl")
using JFinEALE.FESetModule
export FESet, FESet0Manifold, FESet1Manifold, FESet2Manifold, FESet3Manifold
export manifdim
export nfense
export count
export getconn!
export setotherdimension!
export setlabel!
export subset
export cat
export updateconn!
export FESetP1
export FESetL2
export FESetT3
export FESetQ4
export FESetH8

include("FENodeSetModule.jl")
using JFinEALE.FENodeSetModule
export FENodeSet
export spacedim
export xyz3
export count

include("NodalFieldModule.jl")
using JFinEALE.NodalFieldModule
export NodalField
export ndofn
export nfens
export gathersysvec
export gathervalues!
export gatherdofnums!
export numberdofs!
export setebc!
export applyebc!
export scattersysvec!
export clone

include("MeshUtilModule.jl")
using JFinEALE.MeshUtilModule

include("MeshSelectionModule.jl")
using JFinEALE.MeshSelectionModule
export connectedfes
export connectednodes
export fenodeselect
export feselect

include("MeshExportModule.jl")
using JFinEALE.MeshExportModule
export vtkexportmesh
export finealemesh
export graphcsv

include("MeshModificationModule.jl")
using JFinEALE.MeshModificationModule
export meshboundary
export fusenodes
export compactfens
export mergemeshes
export mergenmeshes!
export mergenodes
export renumberfeconn

include("MeshTransformModule.jl")
using JFinEALE.MeshTransformModule
export mirrormesh

include("MeshQuadrilateralModule.jl")
using JFinEALE.MeshQuadrilateralModule
export Q4annulus
export Q4quadrilateral
export Q4elliphole
export Q4block
export Q4blockx
include("MeshLineModule.jl")
using JFinEALE.MeshLineModule
export L2block
export L2blockx

include("MeshTriangleModule.jl")
using JFinEALE.MeshTriangleModule
export T3block

include("MeshHexahedronModule.jl")
using JFinEALE.MeshHexahedronModule
export H8block
export H8blockx
export H8sphere
export H8refine
export H8toH27
export H8hexahedron
export H27block
export H8extrudeQ4
export H8spheren

include("ForceIntensityModule.jl")
using JFinEALE.ForceIntensityModule
export ForceIntensity
export getforce!

include("PropertyHeatDiffusionModule.jl")
using JFinEALE.PropertyHeatDiffusionModule
export PropertyHeatDiffusion

include("MaterialHeatDiffusionModule.jl")
using JFinEALE.MaterialHeatDiffusionModule
export MaterialHeatDiffusion

include("PropertyAcousticFluidModule.jl")
using JFinEALE.PropertyAcousticFluidModule
export PropertyAcousticFluid

include("MaterialAcousticFluidModule.jl")
export MaterialAcousticFluid

using JFinEALE.MaterialAcousticFluidModule
export MaterialAcousticFluid

include("AssemblyModule.jl")
using JFinEALE.AssemblyModule
export SysmatAssemblerBase, SysmatAssemblerSparse, startassembly!, assemble!, makematrix!, SysmatAssemblerSparseSymm, SysvecAssemblerBase, makevector!

include("IntegRuleModule.jl")
using JFinEALE.IntegRuleModule
export IntegRule
export TriRule
export GaussRule

include("FEMMBaseModule.jl")
using JFinEALE.FEMMBaseModule
export FEMMBase
export setRm!
export getRm!
export integrationdata
export integratefieldfunction
export integratefunction
export distribloads
export distribloads

include("FEMMHeatDiffusionModule.jl")
using JFinEALE.FEMMHeatDiffusionModule
export FEMMHeatDiffusion
export conductivity
export conductivity
export nzebcloadsconductivity
export nzebcloadsconductivity

include("FEMMAcousticsModule.jl")
using JFinEALE.FEMMAcousticsModule
export FEMMAcoustics, acousticmass, nzebcloadsacousticmass, acousticstiffness, nzebcloadsacousticstiffness, acousticABC

include("AcousticsAlgorithmModule.jl")
using JFinEALE.AcousticsAlgorithmModule

include("HeatDiffusionAlgorithmModule.jl")
using JFinEALE.HeatDiffusionAlgorithmModule

end


