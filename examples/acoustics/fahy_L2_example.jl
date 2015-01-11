using FESetModule
using MeshLineModule
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

1D mesh.
""")

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
graphics= true;# plot the solution as it is computed?
n=40;#
neigvs=8;
OmegaShift=10.0;

fens,fes = L2block(L,n); # Mesh
setotherdimension!(fes,A)


geom = NodalField(name ="geom",data =fens.xyz)
P = NodalField(name ="P",data =zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoustics(FEMMBase(fes, GaussRule(order=2,dim=1)),
                     MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho)))

S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv =eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")


println("Total time elapsed = ",time() - t0,"s")

using Winston
en=2
pl = FramedPlot(title="Fahy example, mode $en",xlabel="x",ylabel="P")
setattr(pl.frame, draw_grid=true)
ix=sortperm(geom.values[:])
add(pl, Curve(geom.values[:][ix],v[:,en][ix], color="blue"))

display(pl)

true
