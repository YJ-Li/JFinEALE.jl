using FESetModule
using FENodeSetModule
using MeshQuadrilateralModule
using MeshSelectionModule
using MeshModificationModule
using NodalFieldModule
using IntegRuleModule
using PropertyHeatDiffusionModule
using MaterialHeatDiffusionModule
using FEMMBaseModule
using FEMMHeatDiffusionModule
using ForceIntensityModule

println("""

Annular region, ingoing and outgoing flux. Minimum/maximum temperature ~(+/-)0.591.
Mesh of linear quadrilaterals.
""")

t0 = time()

kappa=0.2*[1.0 0; 0 1.0]; # conductivity matrix
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

geom = NodalField(name ="geom",data =fens.xyz)
Temp = NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))


l1 =fenodeselect(fens; box=[0.0 0.0 -rex -rex], inflate = tolerance)
setebc!(Temp,l1,[true],[1],[0.0])
applyebc!(Temp)

numberdofs!(Temp)


p=PropertyHeatDiffusion(kappa)
material=MaterialHeatDiffusion (p)
femm = FEMMHeatDiffusion(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

@time K=conductivity(femm, geom, Temp)

l1=feselect(fens,edge_fes,box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
el1femm = FEMMBase(FESetModule.subset(edge_fes,l1), GaussRule(order=2,dim=1))
fi = ForceIntensity(-magn);#entering the domain
@time F1=-1.0* distribloads(el1femm, geom, Temp, fi, 2);

l1=MeshSelectionModule.feselect(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
el1femm =  FEMMBase(FESetModule.subset(edge_fes,l1), GaussRule(order=2,dim=1))
fi = ForceIntensity(+magn);#leaving the domain
@time F2= -1.0* distribloads(el1femm, geom, Temp, fi, 2);

@time F3=nzebcloadsconductivity(femm, geom, Temp);


@time K=cholfact(K)
@time U=  K\(F1+F2+F3)
@time scattersysvec!(Temp,U[:])

println("Total time elapsed = ",time() - t0,"s")

using MeshExportModule

File =  "annulus.vtk"
vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

true
