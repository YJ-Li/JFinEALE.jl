using JFFoundationModule
using FESetModule
using MeshExportModule
using MeshQuadrilateralModule
using IntegRuleModule: GaussRule
using MeshModificationModule
using MeshSelectionModule
using PropertyHeatDiffusionModule
using MaterialHeatDiffusionModule
using FEMMBaseModule
using FEMMHeatDiffusionModule
using HeatDiffusionAlgorithmModule

println("""

Annular region, ingoing and outgoing flux. Temperature at one node prescribed. 
Minimum/maximum temperature ~(+/-)0.591.
Mesh of linear quadrilaterals.
This version uses the JFinEALE algorithm module.
""")

t0 = time()

kappa=0.2*[1.0 0.0; 0.0 1.0]; # conductivity matrix
magn = 0.06;# heat flux along the boundary
rin= 1.0;#internal radius
rex= 2.0;#external radius
nr=10; nc=80;
Angle=2*pi;
thickness= 1.0;
tolerance=min(rin/nr, rin/nc/2/pi)/10000;



fens,fes = Q4annulus(rin,rex,nr,nc,Angle)
setotherdimension!(fes,thickness)
fens,fes = mergenodes(fens, fes, tolerance);
edge_fes = meshboundary (fes);
setotherdimension!(edge_fes,thickness)

# At a single point apply an essential boundary condition (pin down the temperature)
l1 =fenodeselect(fens; box=[0.0 0.0 -rex -rex], inflate = tolerance)
essential1=dmake(node_list=l1,temperature=0.0)

# The flux boundary condition is applied at two pieces of surface
# Side 1
l1=feselect(fens,edge_fes,box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
el1femm = FEMMBase(subset(edge_fes,l1),GaussRule(order=2,dim=1))
flux1=dmake(femm=el1femm,normal_flux=-magn) # entering the domain
# Side 2
l2=feselect(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
el2femm = FEMMBase(subset(edge_fes,l2),GaussRule(order=2,dim=1))
flux2=dmake(femm=el2femm,normal_flux=+magn) # leaving the domain
 
material=MaterialHeatDiffusion (PropertyHeatDiffusion(kappa))
femm = FEMMHeatDiffusion(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

# Make model data
modeldata= dmake(fens= fens,
                 region=[dmake(femm=femm)],
                 boundary_conditions=dmake(flux=[flux1,flux2],essential=[essential1]));

# Call the solver
modeldata=HeatDiffusionAlgorithmModule.steadystate(modeldata)
geom=modeldata["geom"]
Temp=modeldata["temp"]
println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

println("Total time elapsed = ",time() - t0,"s")

# Postprocessing
MeshExportModule.vtkexportmesh ("annulusmod.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")
