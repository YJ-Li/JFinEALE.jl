using FESetModule
using MeshSelectionModule
using MeshModificationModule
using MeshLineModule
using MeshHexahedronModule
using NodalFieldModule
using IntegRuleModule
using PropertyAcousticFluidModule
using MaterialAcousticFluidModule
using FEMMBaseModule
using FEMMAcousticsModule
using ForceIntensityModule
using PhysicalUnitModule
phun=PhysicalUnitModule.phun
using Base.Test

function fahyL2example()


  println("""

  Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
  Frank J. Fahy, Paolo Gardonio, page 483.

  1D mesh of L2 elements.
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

  @test (fs[2]-343.08)/343.08 < 0.001

  println("Total time elapsed = ",time() - t0,"s")

  # using Winston
  # en=2
  # pl = FramedPlot(title="Fahy example, mode $en",xlabel="x",ylabel="P")
  # setattr(pl.frame, draw_grid=true)
  # ix=sortperm(geom.values[:])
  # add(pl, Curve(geom.values[:][ix],v[:,en][ix], color="blue"))

  # display(pl)

  true
end
fahyL2example()

function fahyH8example()
  println("""

  Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
  Frank J. Fahy, Paolo Gardonio, page 483.

  Hexahedral mesh of H8 elements.
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


  @test (fs[2]-343.08)/343.08 < 0.001

  # File =  "fahy_H8.vtk"
  # en=5;
  # vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars=v[:,en], scalars_name ="Pressure_mode_$n")

  true
end
fahyH8example()

function fahyH27example()
  println("""
  
  Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
  Frank J. Fahy, Paolo Gardonio, page 483.

  Hexahedral mesh of quadratic H27 elements.
  """)

  t0 = time()

  rho=1.21*1e-9;# mass density
  c =343.0*1000;# millimeters per second
  bulk= c^2*rho;
  L=500.0;# length of the box, millimeters
  A=200.0; # cross-sectional area of the box
  graphics= true;# plot the solution as it is computed?
  n=4;#
  neigvs=8;
  OmegaShift=1000.0;

  fens,fes = H8block(L,sqrt(A),sqrt(A),n,1,1); # Mesh
  fens,fes = H8toH27(fens,fes)
  # Alternative meshing
  #fens,fes = H8hexahedron([0.0 0.0 0.0; L sqrt(A) sqrt(A)],n,1,1;block_mesh_handle=H27block)

  geom = NodalField(name ="geom",data =fens.xyz)
  P = NodalField(name ="P",data =zeros(size(fens.xyz,1),1))

  numberdofs!(P)

  femm = FEMMAcoustics(FEMMBase(fes, GaussRule(order=3,dim=3)),
  MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho)))

  S = acousticstiffness(femm, geom, P);
  C = acousticmass(femm, geom, P);

  d,v,nev,nconv =eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
  d = d - OmegaShift;
  fs=real(sqrt(complex(d)))/(2*pi)
  println("Eigenvalues: $fs [Hz]")


  @test (fs[2]-343.08)/343.08 < 0.001

  println("Total time elapsed = ",time() - t0,"s")

  # using MeshExportModule

  # File =  "fahy_H8.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

  # println(" Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

  true
end
fahyH27example()
