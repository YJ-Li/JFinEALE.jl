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

include("MaterialOrientationModule.jl") 
using JFinEALE.MaterialOrientationModule
export MaterialOrientation
export updateRm!

include("FESetModule.jl")
using JFinEALE.FESetModule
export FESet, FESet0Manifold, FESet1Manifold, FESet2Manifold, FESet3Manifold
export manifdim
export nfensperfe
export count
export getlabel
export getconn!
export setotherdimension!
export setlabel!
export subset
export cat
export updateconn!
export FESetP1
export FESetL2
export FESetL3
export FESetT3
export FESetQ4
export FESetQ9
export FESetQ8
export FESetT6
export FESetH8
export FESetH20
export FESetH27
export FESetT4
export FESetT10


include("FENodeSetModule.jl")
using JFinEALE.FENodeSetModule
export FENodeSet
export spacedim
export xyz3
export count

include("FENodeToFEMapModule.jl")
using JFinEALE.FENodeToFEMapModule
export FENodeToFEMap

include("NodalFieldModule.jl")
using JFinEALE.NodalFieldModule
export NodalField
export ndofn
export nfens
export gathersysvec
export gathervaluesasvec!
export gathervaluesasmat!
export gatherdofnumsasvec!
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
export Q8block
export Q4toQ8


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

include("MeshTetrahedronModule.jl")
using JFinEALE.MeshTetrahedronModule
export T4blocka
export T4blockb
export T4blockca
export T4blockcb
export T4block
export T4blockx
export T4toT10
export T10block

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
using JFinEALE.MaterialAcousticFluidModule
export MaterialAcousticFluid

include("PropertyDeformationLinearModule.jl")
using JFinEALE.PropertyDeformationLinearModule
export PropertyDeformationLinear
export tangentmoduli3d!

include("PropertyDeformationLinearIsoModule.jl")
using JFinEALE.PropertyDeformationLinearIsoModule
export PropertyDeformationLinearIso
export tangentmoduli3d!

include("PropertyDeformationLinearOrthoModule.jl")
using JFinEALE.PropertyDeformationLinearOrthoModule
export PropertyDeformationLinearOrtho
export tangentmoduli3d!

include("DeformationModelReductionModule.jl")
using JFinEALE.DeformationModelReductionModule
export DeformationModelReduction
export DeformationModelReduction1D
export DeformationModelReduction2DStrain
export DeformationModelReduction2DStress
export DeformationModelReduction2DAxisymm
export DeformationModelReduction3D
export nstrains
export Blmat!

include("MaterialDeformationModule.jl")
using JFinEALE.MaterialDeformationModule
export strain2x2tto3v!
export strain3vto2x2t!
export strain3x3tto6v!
export strain6vto3x3t!
export stress2x2to3v!
export stress3vto2x2t!
export stress3vto3x3t!
export stress4vto3x3t!
export stress6vto3x3t!
export stress3x3tto6v!
export strain9vto6v!
export strain6vto9v!
export stress9vto6v!
export stressvectorrotation
export strainvectorrotation
export stressvectorrotation
export strainvectorrotation
export rotatestiffness!
export rotatecompliance!


include("MaterialDeformationLinearModule.jl")
using JFinEALE.MaterialDeformationLinearModule
export MaterialDeformationLinear
export tangentmoduli!
export update!
export thermalstress

include("AssemblyModule.jl")
using JFinEALE.AssemblyModule
export SysmatAssemblerBase
export SysmatAssemblerSparse
export startassembly!
export assemble!
export makematrix!
export SysmatAssemblerSparseSymm
export startassembly!
export assemble!
export makematrix!
export SysvecAssemblerBase
export SysvecAssembler
export startassembly!
export assemble!
export makevector!


include("IntegRuleModule.jl")
using JFinEALE.IntegRuleModule
export IntegRule
export TriRule
export GaussRule
export TetRule

include("FEMMBaseModule.jl")
using JFinEALE.FEMMBaseModule
export FEMMBase
export setRm!
export getRm!
export integrationdata
export integratefieldfunction
export integratefunction
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
export FEMMAcoustics
export acousticmass
export nzebcloadsacousticmass
export acousticstiffness
export nzebcloadsacousticstiffness
export acousticABC


include("FEMMDeformationLinearModule.jl")
using JFinEALE.FEMMDeformationLinearModule
export FEMMDeformationLinear
export stiffness
export nzebcloadsstiffness
export mass
export inspectintegpoints
export fieldfromintegpoints


include("AcousticsAlgorithmModule.jl")
using JFinEALE.AcousticsAlgorithmModule

include("HeatDiffusionAlgorithmModule.jl")
using JFinEALE.HeatDiffusionAlgorithmModule

end


