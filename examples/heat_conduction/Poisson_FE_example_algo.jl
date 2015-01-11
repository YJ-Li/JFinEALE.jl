using JFFoundationModule
using MeshExportModule
using IntegRuleModule
using MeshSelectionModule
using MeshTriangleModule
using HeatDiffusionAlgorithmModule
using FEMMBaseModule




A= 1.0
thermal_conductivity=eye(2,2); # conductivity matrix
magn = -6.0; #heat source
boundaryf(x)=1.0 + x[:,1].^2 + 2.0*x[:,2].^2;
N=20;

println("""

Heat conduction example described by Amuthan A. Ramabathiran
http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
Unit square, with known temperature distribution along the boundary, 
and uniform heat generation rate inside.  Mesh of regular TRIANGLES,
in a grid of $N x $N edges. 
This version uses the JFinEALE algorithm module.
"""
)
t0 = time()

fens,fes =T3block(A, A, N, N)


# Define boundary conditions
l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)

essential1= dmake(node_list=[l1 l2 l3 l4],temperature=boundaryf);

# Make model data
modeldata= dmake(fens= fens,
                 region=[dmake(conductivity=thermal_conductivity,
                               Q=(x,J)->[magn],fes=fes,integration_rule=TriRule(npts=1))],
                 boundary_conditions=dmake(essential=[essential1]));


# Call the solver
modeldata=HeatDiffusionAlgorithmModule.steadystate(modeldata)

println("Total time elapsed = ",time() - t0,"s")

geom=modeldata["geom"]
Temp=modeldata["temp"]
femm=modeldata["region"][1]["femm"]
function errfh(loc,val)
    exact = boundaryf(loc)
    return ((exact-val)*exact)[1]
end

femm.femmbase.integration_rule=TriRule(npts=6)
E= integratefieldfunction (femm.femmbase,geom,Temp,errfh,0.0,m=3)
println("Error=$E")                                             

# Postprocessing
# geom=modeldata["geom"]
# Temp=modeldata["temp"]
# MeshExportModule.vtkexportmesh ("a.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

