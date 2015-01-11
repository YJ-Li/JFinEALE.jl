using FESetModule
using FENodeSetModule
using MeshTriangleModule
using MeshSelectionModule
using MeshModificationModule
using NodalFieldModule
using IntegRuleModule
using PropertyAcousticFluidModule
using MaterialAcousticFluidModule
using FEMMBaseModule
using FEMMAcousticsModule
using ForceIntensityModule
using PhysicalUnitModule
phun=PhysicalUnitModule.phun

t0 = time()

rho=1.21*phun("kg/m^3");# mass density
c =343.0*phun("m/s");# sound speed
bulk= c^2*rho;
omega= 54.5901*phun("rev/s")
vn0= -1.0*phun("m/s")
Lx=10.0*phun("m");# length of the box, millimeters
Ly=1.0*phun("m"); # length of the box, millimeters
n=100;#number of elements along the length

println("""

Straight duct with anechoic termination.
Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44. 
Both real and imaginary components of the pressure should have amplitude of
rho*c=$(rho*c).

Triangle mesh. 
""")

fens,fes = T3block(Lx,Ly,n,2); # Mesh
bfes = meshboundary(fes)
L0=feselect(fens,bfes,facing=true,direction=[-1.0 0.0])
L10=feselect(fens,bfes,facing=true,direction=[+1.0 0.0])
Lx=feselect(fens,bfes,facing=true,direction=[0.0 +1.0])
nLx=connectednodes(FESetModule.subset(bfes,Lx))

geom = NodalField(name ="geom",data =fens.xyz)
P = NodalField(name ="P",data =zeros(Complex128,size(fens.xyz,1),1))

numberdofs!(P)

material=MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho))
femm = FEMMAcoustics(FEMMBase(fes, TriRule(npts=1)), material)
 
S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);


E10femm = FEMMAcoustics(FEMMBase(subset(bfes,L10), GaussRule(order=2,dim=1)), material)
D = acousticABC(E10femm, geom, P);

E0femm = FEMMBase(subset(bfes,L0), GaussRule(order=2,dim=1))
fi = ForceIntensity(-1.0im*omega*rho*vn0);
F = distribloads(E0femm, geom, P, fi, 2);

p=(-omega^2*S +omega*1.0im*D + C)\F              # 
scattersysvec!(P,p[:])


println("Pressure amplitude bounds: ")
println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")

println("Total time elapsed = ",time() - t0,"s")

using MeshExportModule

# File =  "straight_duct.vtk"
# scalars=real(P.values);
# MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values scalars/100.0], MeshExportModule.T3; scalars=scalars, scalars_name ="Pressure")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)


# ix=sortperm(geom.values[nLx,1])

# MeshExportModule.graphcsv(geom.values[nLx,1][ix],scalars[nLx][ix];xlabel="x",ylabel="Pressure",file= "agraph")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" agraph.csv`)

scalars=real(P.values);
ix=sortperm(geom.values[nLx,1])
# using Gadfly                    # 
# pl=plot(x=geom.values[nLx,1][ix],y=scalars[nLx][ix], Geom.line)
# draw(PNG("myplot.png", 12cm, 6cm), pl)

using Winston

pl = FramedPlot(title="Straight duct with anechoic termination",xlabel="x",ylabel="Re P, Im P")
setattr(pl.frame, draw_grid=true)
ix=sortperm(geom.values[nLx,1])
add(pl, Curve(geom.values[nLx,1][ix],real(P.values)[nLx][ix], color="blue"))
add(pl, Curve(geom.values[nLx,1][ix],imag(P.values)[nLx][ix], color="red"))
add(pl, Curve(geom.values[:,1],imag(P.values), color="green"))

# pl=plot(geom.values[nLx,1][ix],scalars[nLx][ix])
# xlabel("x")
# ylabel("Pressure")
display(pl)


true
