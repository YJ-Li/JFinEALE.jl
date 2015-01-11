using FESetModule
using MeshQuadrilateralModule
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

println("""

Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 123. 
Internal resonance problem. Reference frequencies: 90.7895, 181.579, 215.625, 233.959, 272.368, 281.895 
Quadrilateral mesh. 
""")

t0 = time()

rho=1.21*1e-9;# mass density
c =345.0*1000;# millimeters per second
bulk= c^2*rho;
Lx=1900.0;# length of the box, millimeters
Ly=800.0; # length of the box, millimeters
n=140;#
neigvs=18;
OmegaShift=10.0;

fens,fes = Q4block(Lx,Ly,n,n); # Mesh

geom = NodalField(name ="geom",data =fens.xyz)
P = NodalField(name ="P",data =zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoustics(FEMMBase(fes, GaussRule(order=2,dim=2)),
                     MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho)))
 
S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv =eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")


println("Total time elapsed = ",time() - t0,"s")

using MeshExportModule

File =  "rigid_box.vtk"
n=15;
scalars=v[:,n];
vtkexportmesh (File, fes.conn, [geom.values scalars], MeshExportModule.Q4; scalars=scalars, scalars_name ="Pressure_mode_$n")
@async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

true
