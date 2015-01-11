using MeshHexahedronModule
using MeshSelectionModule
using MeshModificationModule
using NodalFieldModule
using IntegRuleModule
using PropertyAcousticFluidModule
using MaterialAcousticFluidModule
using FEMMAcousticsModule
using ForceIntensityModule
using PhysicalUnitModule
phun=PhysicalUnitModule.phun

println("""
Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
Frank J. Fahy, Paolo Gardonio, page 483. 

Hexahedral mesh.
""")

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
n=40;#
neigvs=8;
OmegaShift=10.0;

fens,fes = H8block(L,sqrt(A),sqrt(A),n,1,1); # Mesh

geom = NodalField(name ="geom",data =fens.xyz)
P = NodalField(name ="P",data =zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoustics(FEMMBase(fes, GaussRule(order=2,dim=3)),
                     MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho)))
                     
 
S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv =eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")


println("Total time elapsed = ",time() - t0,"s")

using MeshExportModule

File =  "fahy_H8.vtk"
en=5;
vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars=v[:,en], scalars_name ="Pressure_mode_$n")

true
