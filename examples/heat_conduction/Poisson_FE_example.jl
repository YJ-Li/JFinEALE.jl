using FESetModule
using FENodeSetModule
using MeshTriangleModule
using MeshSelectionModule
using NodalFieldModule; NF=NodalFieldModule;
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
and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom). 
"""
)
t0 = time()

A= 1.0 # dimension of the domain (length of the side of the square)
thermal_conductivity= eye(2,2); # conductivity matrix
magn = -6.0; #heat source
tempf(x)=(1.0 + x[:,1].^2 + 2*x[:,2].^2);#the exact distribution of temperature
N=1000;# number of subdivisions along the sides of the square domain


println("Mesh generation")
@time fens,fes =T3block(A, A, N, N)

geom = NF.NodalField(name ="geom",data =fens.xyz)
Temp = NF.NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))

println("Searching nodes  for BC")
@time l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
@time l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
@time l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
@time l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)
List=[l1, l2, l3, l4];
NF.setebc!(Temp,List,trues(length(List)),List*0+1,tempf(geom.values[List,:])[:])
NF.applyebc!(Temp)
NF.numberdofs!(Temp)

t1 = time()

p=PropertyHeatDiffusion(thermal_conductivity)
material=MaterialHeatDiffusion (p)

femm = FEMMHeatDiffusion(FEMMBase(fes, TriRule(npts=1)), material)


println("Conductivity")
@time K=conductivity(femm, geom, Temp)
println("Nonzero EBC")
@time F2=nzebcloadsconductivity(femm, geom, Temp);
println("Internal heat generation")
fi = ForceIntensity(magn);
@time F1=distribloads(femm.femmbase, geom, Temp, fi, 3);

println("Factorization")
@time K=cholfact(K)
println("Solution of the factorized system")
@time U=  K\(F1+F2)
NF.scattersysvec!(Temp,U[:])

println("Total time elapsed = $(time() - t0) [s]")
println("Solution time elapsed = $(time() - t1) [s]")

Error= 0.0
for k=1:size(fens.xyz,1)
    Error=Error+abs(Temp.values[k,1]-tempf(fens.xyz[k,:]))
end
println("Error =$Error")

# using MeshExportModule

# File =  "a.vtk"
# MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

true
