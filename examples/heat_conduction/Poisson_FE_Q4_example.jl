using FESetModule
using FENodeSetModule
using MeshQuadrilateralModule
using MeshSelectionModule
using NodalFieldModule
using IntegRuleModule
using PropertyHeatDiffusionModule
using MaterialHeatDiffusionModule
using FEMMBaseModule
using FEMMHeatDiffusionModule
using ForceIntensityModule

println("""

Heat conduction example described by Amuthan A. Ramabathiran
http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
Unit square, with known temperature distribution along the boundary, 
and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom). 
"""
)
t0 = time()

A= 1.0
thermal_conductivity= eye(2,2); # conductivity matrix
magn = -6.0; #heat source
tempf(x)=(1.0 + x[:,1].^2 + 2*x[:,2].^2);
N=1000;

println("Mesh generation")
@time fens,fes =Q4block(A, A, N, N)



geom = NodalField(name ="geom",data =fens.xyz)
Temp = NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))


println("Searching nodes  for BC")
@time l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
@time l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
@time l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
@time l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)
List=[l1, l2, l3, l4];
setebc!(Temp,List,trues(length(List)),List*0+1,tempf(geom.values[List,:])[:])
applyebc!(Temp)

numberdofs!(Temp)

t1 = time()
 
m=MaterialHeatDiffusion (PropertyHeatDiffusion(thermal_conductivity))
femm = FEMMHeatDiffusion(FEMMBase(fes, GaussRule(order=2,dim=2)), m)

println("Conductivity")
@time K=conductivity(femm, geom, Temp)
#Profile.print()

println("Nonzero EBC")
@time F2=nzebcloadsconductivity(femm, geom, Temp);
println("Internal heat generation")
fi = ForceIntensity(magn);
@time F1=distribloads(femm.femmbase, geom, Temp, fi, 3);


println("Factorization")
@time K=cholfact(K)
println("Solution of the factorized system")
@time U=  K\(F1+F2)
scattersysvec!(Temp,U[:])


println("Total time elapsed = $(time() - t0) [s]")
println("Solution time elapsed = $(time() - t1) [s]")

# using MeshExportModule

# File =  "a.vtk"
# MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

Error= 0.0
for k=1:size(fens.xyz,1)
    Error=Error+abs(Temp.values[k,1]-tempf(fens.xyz[k,:]))
end
println("Error =$Error")


true
