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


    # #println("""

    #         Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    #         Frank J. Fahy, Paolo Gardonio, page 483.

    #         1D mesh of L2 elements.
    #         """)

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
    #println("Eigenvalues: $fs [Hz]")

    @test (fs[2]-343.08)/343.08 < 0.001

    #println("Total time elapsed = ",time() - t0,"s")

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
    # #println("""

    #         Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    #         Frank J. Fahy, Paolo Gardonio, page 483.

    #         Hexahedral mesh of H8 elements.
    #         """)

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
    #println("Eigenvalues: $fs [Hz]")


    #println("Total time elapsed = ",time() - t0,"s")


    @test (fs[2]-343.08)/343.08 < 0.001

    # File =  "fahy_H8.vtk"
    # en=5;
    # vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars=v[:,en], scalars_name ="Pressure_mode_$n")

    true
end
fahyH8example()

function fahyH27example()
    # println("""
            
    #         Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    #         Frank J. Fahy, Paolo Gardonio, page 483.

    #         Hexahedral mesh of quadratic H27 elements.
    #         """)

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
    #println("Eigenvalues: $fs [Hz]")


    @test (fs[2]-343.08)/343.08 < 0.001

    #println("Total time elapsed = ",time() - t0,"s")

    # using MeshExportModule

    # File =  "fahy_H8.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

    # #println(" Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

    true
end
fahyH27example()

using FESetModule
using FENodeSetModule
using MeshHexahedronModule
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
function straight_duct_H8_1()



    t0 = time()

    rho=1.21*phun("kg/m^3");# mass density
    c =343.0*phun("m/s");# sound speed
    bulk= c^2*rho;
    omega= 54.5901*phun("rev/s")
    vn0= -1.0*phun("m/s")
    Lx=10.0*phun("m");# length of the box, millimeters
    Ly=1.0*phun("m"); # length of the box, millimeters
    n=20;#number of elements along the length 

    # println("""

    #         Straight duct with anechoic termination.
    #         Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44. 
    #         Both real and imaginary components of the pressure should have amplitude of
    #         rho*c=$(rho*c).

    #         Hexahedral mesh. 
    #         """)

    fens,fes = H8block(Lx,Ly,Ly,n,2,2); # Mesh
    bfes = meshboundary(fes)
    L0=feselect(fens,bfes,facing=true,direction=[-1.0 0.0 0.0])
    L10=feselect(fens,bfes,facing=true,direction=[+1.0 0.0 0.0])
    nLx=fenodeselect(fens,box=[0.0 Lx  0.0 0.0 0.0 0.0],inflate=Lx/1.0e5)

    geom = NodalField(name ="geom",data =fens.xyz)
    P = NodalField(name ="P",data =zeros(Complex128,size(fens.xyz,1),1))

    numberdofs!(P)


    material=MaterialAcousticFluid (PropertyAcousticFluid(bulk,rho))
    femm = FEMMAcoustics(FEMMBase(fes, GaussRule(order=2,dim=3)), material)
    
    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);


    E10femm = FEMMAcoustics(FEMMBase(subset(bfes,L10),GaussRule(order=2,dim=2)), material)
    D = acousticABC(E10femm, geom, P);

    E0femm = FEMMBase(subset(bfes,L0), GaussRule(order=2,dim=2))
    fi = ForceIntensity(-1.0im*omega*rho*vn0);
    F = distribloads(E0femm, geom, P, fi, 2);

    p=(-omega^2*S +omega*1.0im*D + C)\F              # 
    scattersysvec!(P,p[:])

    #println("Pressure amplitude bounds: ")
    #println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    #println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")

    ref=rho*c
    @test abs(minimum(real(P.values))-(-ref))/ref < 0.02
    @test abs(minimum(imag(P.values))-(-ref))/ref < 0.02
    @test abs(maximum(real(P.values))-(+ref))/ref < 0.02
    @test abs(maximum(imag(P.values))-(+ref))/ref < 0.02
    #println("Total time elapsed = ",time() - t0,"s")

    # using MeshExportModule

    # File =  "straight_duct.vtk"
    # scalars=real(P.values);
    # MeshExportModule.vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars=scalars, scalars_name ="Pressure")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)


    # ix=sortperm(geom.values[nLx,1])

    #scalars=real(P.values);
    # MeshExportModule.graphcsv(geom.values[nLx,1][ix],scalars[nLx][ix];xlabel="x",ylabel="Pressure",file= "agraph")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" agraph.csv`)

    # scalars=real(P.values);
    # ix=sortperm(geom.values[nLx,1])
    # using Gadfly                    # 
    # pl=plot(x=geom.values[nLx,1][ix],y=scalars[nLx][ix], Geom.line)
    # draw(PNG("myplot.png", 12cm, 6cm), pl)

    # using Winston

    # pl = FramedPlot(title="Straight duct with anechoic termination",xlabel="x",ylabel="Re P, Im P")
    # setattr(pl.frame, draw_grid=true)
    # ix=sortperm(geom.values[nLx,1])
    # add(pl, Curve(geom.values[nLx,1][ix],real(P.values)[nLx][ix], color="blue"))
    # add(pl, Curve(geom.values[nLx,1][ix],imag(P.values)[nLx][ix], color="red"))
    # add(pl, Curve(geom.values[:,1],imag(P.values), color="green"))

    # # pl=plot(geom.values[nLx,1][ix],scalars[nLx][ix])
    # # xlabel("x")
    # # ylabel("Pressure")
    # display(pl)

    true
end
straight_duct_H8_1()
