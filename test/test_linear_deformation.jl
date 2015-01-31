using JFinEALE
using Base.Test


function cookstress()
    #println("Cook membrane problem, plane stress."        )        
    t0 = time()

    E=1.0;
    nu=1.0/3;
    width =48.0; height = 44.0; thickness  = 1.0;
    free_height  = 16.0;
    Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
    magn=1./free_height;# Magnitude of applied load
    convutip=23.97;
    n=32;#*int(round(sqrt(170.)/2.)); # number of elements per side
    tolerance=minimum([width,height])/n/1000.;#Geometrical tolerance

    fens,fes =T3block(width,height, n, n)
    setotherdimension!(fes,1.0)

    # Reshape into a trapezoidal panel
    for i=1:count(fens)
        fens.xyz[i,2]=fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
    end

    geom = NodalField(name ="geom",data =fens.xyz)
    u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

    l1 =fenodeselect(fens; box=[0,0,-Inf, Inf], inflate = tolerance)
    setebc!(u,l1,trues(length(l1)),l1*0+1,[0.0])
    setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    Toplist  =feselect(fens,boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
    el1femm =  FEMMBase(subset(boundaryfes,Toplist), GaussRule(order=2,dim=1))
    fi = ForceIntensity([0.0,+magn]);
    F2= distribloads(el1femm, geom, u, fi, 2);


    p=PropertyDeformationLinearIso(E,nu)
    material=MaterialDeformationLinear (p)

    femm = FEMMDeformationLinear(FEMMBase(fes, TriRule(npts=1)), material)

    K =stiffness(DeformationModelReduction2DStress, femm, geom, u)
    #K=cholfact(K)
    U=  K\(F2)
    scattersysvec!(u,U[:])

    nl=fenodeselect (fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],inflate=tolerance);
    theutip=zeros(JFFlt,1,2)
    gathervaluesasmat!(u,theutip,nl);
    #    println("$(time()-t0) [s];  displacement =$(theutip[2]) as compared to converged $convutip") 

    @test abs(theutip[2]-23.8155)/23.8155 < 1.e-3 # FinEALE solution
end
cookstress()


module scratch1

using JFinEALE
using Base.Test

## Thick pipe with internal pressure: axially symmetric model
#

##
# Link to the  <matlab:edit('pub_thick_pipe_axi') m-file>.

## Description
##
# This is a simple modification of the full three-dimensional simulation of
# the tutorial pub_thick_pipe that implements the axially-symmetric model
# reduction procedure.
##
# An infinitely long thick walled cylindrical pipe
# with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
# subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
# 2 mm and a 90-degree angle sector is considered for the finite element
# analysis. The material properties are taken as  isotropic linear elastic
# with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
# behavior. This problem has been proposed to by MacNeal and Harder as a
# test of an element's ability to represent the  response of a nearly
# incompressible material. The plane-strain condition is assumed in the
# axial direction of the pipe which together with the radial symmetry
# confines the material in all but the radial direction and therefore
# amplifies the numerical difficulties associated with the confinement of
# the nearly incompressible material.
##
# There is an analytical solution to this problem. Timoshenko and Goodier
# presented the original solution of Lame in their textbook. We are going
# to compare with  both the stress distribution (radial and hoop stresses)
# and the displacement of the inner  cylindrical surface.

##
#
# <html>
# <table border=0><tr><td>
# <img src="../docs/pub_thick_pipe_axi.png" width = "30#">
# </td></tr>
# <tr><td>Figure 1. Definition of the geometry of the internally pressurized thick pipe</td></tr>
# </table>
# </html>

##
# References:
#
# # Macneal RH, Harder RL (1985) A proposed standard set of problems to test
# finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
#
# # Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

## Solution
#

type MyIData
    c::JFInt
    r::JFFltVec
    s::JFFltVec
end

function pub_thick_pipe_axi()

    #using JFinEALE

    ##
    # Internal radius of the pipe.
    a=3*phun("MM");
    ##
    # External radius of the pipe.
    b=9*phun("MM");
    ##
    # Thickness of the slice.
    t=2*phun("MM");

    ##
    # Geometrical tolerance.
    tolerance  =a/10000.;
    ##
    # Young's modulus and Poisson's ratio.
    E=1000*phun("MEGA*PA");
    nu=0.499;
    ##
    # Applied pressure on the internal surface.
    press=  1.0*phun("MEGA*PA");

    ##
    # Analytical solutions.   Radial stress:
    radial_stress(r) =press*a.^2/(b^2-a^2).*(1-b^2./r.^2);
    ##
    # Circumferential (hoop) stress:
    hoop_stress(r)=press*a.^2/(b^2-a^2).*(1+b^2./r.^2);

    ##
    # Radial displacement:
    radial_displacement(r)=press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

    ##
    # Therefore the radial displacement of the loaded surface will be:
    urex = radial_displacement(a);


    ##
    # The mesh parameters: The numbers of element edges axially,
    # and through the thickness of the pipe wall (radially).

    na=1; nt=10;

    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    mr=DeformationModelReduction2DAxisymm


    # Create the mesh and initialize the geometry.  First we are going
    # to construct the block of elements with the first coordinate
    # corresponding to the thickness in the radial direction, and the second
    # coordinate is the thickness in the axial direction.
    fens,fes =  Q8block (b-a, t, nt, na);
    fes.axisymm=(mr== DeformationModelReduction2DAxisymm)           # note that this reflects the chosen model reduction

    # Extract the boundary  and mark the finite elements on the
    # interior surface.
    bdryfes = meshboundary(fes);
    bdryfes.axisymm=(mr== DeformationModelReduction2DAxisymm)

    bcl = feselect(fens, bdryfes, box=[0.,0.,-Inf,Inf],inflate=tolerance);
    internal_fenids= connectednodes (subset(bdryfes,bcl));
    # Now  shape the block  into  the actual wedge piece of the pipe.
    for i=1:count (fens)
        fens.xyz[i,:]=fens.xyz[i,:]+[a 0.0];
    end

    # now we create the geometry and displacement fields
    geom = NodalField(name ="geom",data =fens.xyz)
    u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

    # The plane-strain condition in the axial direction  is specified by selecting nodes
    # on the plane y=0 and y=t.
    l1 =fenodeselect(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
    setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])
    l1 =fenodeselect(fens; box=[-Inf Inf t t], inflate = tolerance)
    setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])    

    applyebc!(u)
    numberdofs!(u)

    # The traction boundary condition is applied in the radial
    # direction.

    el1femm =  FEMMBase(subset(bdryfes,bcl), GaussRule(order=3,dim=1))
    fi = ForceIntensity(reshape([press; 0.0],2,1));
    F2= distribloads(el1femm, geom, u, fi, 2);

    # Property and material
    material=MaterialDeformationLinear (PropertyDeformationLinearIso(E,nu))

    femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

    K =stiffness(mr, femm, geom, u)
    #K=cholfact(K)
    U=  K\(F2)
    scattersysvec!(u,U[:])

    # Transfer the solution of the displacement to the nodes on the
    # internal cylindrical surface and convert to
    # cylindrical-coordinate displacements there.
    uv=u.values[internal_fenids,:]
    # Report the  relative displacement on the internal surface:
    #println("(Approximate/true displacement) at the internal surface: $( mean(uv[:,1])/urex*100  ) %")
    @test abs(mean(uv[:,1])/urex-1.0)<1.0e-9
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(mr, femm, geom, u, :Cauchy, 1)


    # File =  "thick_pipe_sigmax.vtk"
    # vtkexportmesh(File, fens, fes; scalars=fld.values,scalars_name ="sigmax")

    # Produce a plot of the solution components in the cylindrical
    # coordinate system.
    # Plot the analytical solution.

    # r  =linspace(a,b,100);
    # using Winston                   # 
    # pl = FramedPlot(title="Thick pipe, axially symmetric solution",xlabel="r",ylabel="Radial stress")
    # setattr(pl.frame, draw_grid=true)
    # add(pl, Curve(r,radial_stress(r), color="black"))
    # display(pl)

    function inspector(idat::MyIData,out,loc,pc)
        push!(idat.r,loc[1])
        push!(idat.s,out[idat.c])
        return idat
    end

    idat=MyIData(1,JFInt[],JFInt[])
    idat=inspectintegpoints(mr, femm, geom, u, [1:count(fes)], inspector, idat; output=:Cauchy)
    #
    #println("$( idat  )")
    
    @test  abs(idat.r[1]-0.003126794919243112)<1.0e-9
    @test  abs(idat.s[1]- -910911.9777008593)<1.0e-2

    ## Discussion
    #
    ##
    # The axially symmetric model is clearly very effective
    # computationally, as the size is much reduced compared to the 3-D
    # model.  In conjunction with  uniform or selective reduced integration
    # it can be very accurate as well.
end
end
scratch1.pub_thick_pipe_axi()

module scratch2
using JFinEALE
using Base.Test


# println("""
#         % Vibration modes of unit cube  of almost incompressible material.
#         %
#         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
#         % tetrahedral. International Journal for Numerical Methods in
#         % Engineering 67: 841-867.""") 
#         t0 = time()


E = 1*phun("PA");
nu = 0.499;
rho= 1*phun("KG/M^3");
a=1*phun("M"); b=a; h= a;
n1=10;# How many element edges per side?
na= n1; nb= n1; nh =n1;
neigvs=20                   # how many eigenvalues
OmegaShift=(0.01*2*pi)^2;

Modelreduction=DeformationModelReduction3D
fens,fes =H20block(a,b,h, na,nb,nh)

geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(size(fens.xyz,1),3)) # displacement field

numberdofs!(u)

p=PropertyDeformationLinearIso(rho,E,nu)
material=MaterialDeformationLinear (p)

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=3)), material)

K =stiffness(Modelreduction, femm, geom, u)
femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=3,dim=3)), material)
M =mass(Modelreduction, femm, geom, u)
d,v,nev,nconv =eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt(complex(d)))/(2*pi)
#println("Eigenvalues: $fs [Hz]")

@test fs[7]-0.262598691962592<1.0e-6

end

module mxxxx1

using JFinEALE
using JFinEALE.DeformationLinearAlgorithmModule
using JFinEALE.MeshExportModule
using Base.Test


# println("""
# % The initially twisted cantilever beam is one of the standard test
# % problems for verifying the finite-element accuracy [1]. The beam is
# % clamped at one end and loaded either with unit in-plane or 
# % unit out-of-plane force at the other. The centroidal axis of the beam is
# % straight at the undeformed  configuration, while its cross-sections are
# % twisted about the centroidal axis from 0 at the clamped end to pi/2 at
# % the free end. 
# %
# % Reference:
# % Zupan D, Saje M (2004) On "A proposed standard set of problems to test
# % finite element accuracy": the twisted beam. Finite Elements in Analysis
# % and Design 40: 1445-1451.  
# """)
function  Twisted_beam(dir)
    E=0.29e8;
    nu=0.22;
    W=1.1;
    L=12.;
    t= 0.32;
    nl=2; nt=1; nw=1; ref=7;
    p=  1/W/t;
    if (dir==2)
        #   Loading in the Y direction
        loadv=[0;p;0];dir=2;uzex=0.001753248285256; # Harder: 1.754e-3;
    else
        #  Loading in the Z direction
        loadv=[0;0;p];dir=3;uzex=0.005424534868469; # Harder: 5.424e-3;
    end
    
    tolerance =t/1000;
    
    fens,fes =H20block(L,W,t, nl*ref,nw*ref,nt*ref)

    # Reshape into a twisted beam shape
    for i=1:count (fens)
        a=fens.xyz[i,1]/L*(pi/2); y=fens.xyz[i,2]-(W/2); z=fens.xyz[i,3]-(t/2);
        fens.xyz[i,:]=[fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    
    # Clamped end of the beam
    l1 =fenodeselect(fens; box=[0 0 -100*W 100*W -100*W 100*W], inflate = tolerance)
    essential1=dmake(node_list=l1,component=1:3,displacement=0.0)

    # Traction on the opposite edge
    boundaryfes =  meshboundary(fes);
    Toplist  =feselect(fens,boundaryfes, box= [L L -100*W 100*W -100*W 100*W], inflate=  tolerance);
    el1femm =  FEMMBase(subset(boundaryfes,Toplist), GaussRule(order=2,dim=2))
    flux1=dmake(femm=el1femm,traction_vector=loadv)

    
    # Make the region
    p=PropertyDeformationLinearIso(E,nu)
    region1=dmake(fes=fes, property=p, integration_rule=GaussRule(order=2,dim=3))

    # Model reduction type
    mr=DeformationModelReduction3D

    # Make model data
    modeldata= dmake(modelreduction=mr,
                     fens= fens,region=[region1],
                     boundary_conditions=dmake(traction=[flux1],essential=[essential1]));

    # Call the solver
    modeldata=JFinEALE.DeformationLinearAlgorithmModule.linearstatics(modeldata)
    geom=modeldata["geom"]
    u=modeldata["u"]

    # Extract the solution
    nl=fenodeselect (fens, box=[L L -100*W 100*W -100*W 100*W],inflate=tolerance);
    theutip=mean(u.values[nl,:],1)
    #println("displacement =$(theutip[dir]) as compared to converged $uzex")
    @test    abs(theutip[dir]-uzex)<1.0e-5
    
    #println("$( abs(theutip[dir]-uzex)  )")
    # # Write out mesh with displacements
    # dadd!(modeldata,postprocessing=dmake(file=  "twisted_beam"))
    # modeldata=JFinEALE.DeformationLinearAlgorithmModule.exportdeformation(modeldata)

    # # Write out mesh with stresses
    # dadd!(modeldata,postprocessing=dmake(file=  "twisted_beam", output=:Cauchy, component=:xy))
    # modeldata=JFinEALE.DeformationLinearAlgorithmModule.exportstress(modeldata)

    # # Write out mesh with von Mises stresses
    # dadd!(modeldata,postprocessing=dmake(file=  "twisted_beam", output=:vm))
    # modeldata=JFinEALE.DeformationLinearAlgorithmModule.exportstress(modeldata)

    true

end
Twisted_beam(2)
Twisted_beam(3)

end

module m # pub_thick_pipe_axi()

using JFinEALE
using Base.Test

# NAFEMS LE11 benchmark with Q4 elements.
# # This is a test recommended by the National Agency for Finite Element
# # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
# # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.  
# #
# # Target solution: Direct stress,  = –105 MPa at point A.
#function  LE11NAFEMS()
# Parameters:
Ea= 210000*phun("MEGA*Pa")
nua= 0.3;
alphaa=2.3e-4;              # thermal expansion coefficient
sigmaA=-105*phun("MEGA*Pa")
nref= 1;                        # how many times should we refine the mesh?
X=[1.     0.;#A
   1.4    0.;#B
   0.995184726672197   0.098017140329561;
   1.393258617341076 0.137223996461385;
   0.980785  0.195090;#
   1.37309939 0.27312645;
   0.956940335732209   0.290284677254462
   1.339716470025092 0.406398548156247
   0.9238795  0.38268;#C
   1.2124  0.7;#D
   0.7071  0.7071;#E
   1.1062  1.045;#F
   0.7071  (0.7071+1.79)/2;#(E+H)/2
   1.      1.39;#G
   0.7071  1.79;#H
   1.      1.79;#I
   ]*phun("M")
   tolerance =1.e-6*phun("M")
   ##
   # Note that the material object needs to be created with the proper
   # model-dimension reduction in mind.  In this case that is the axial symmetry
   # assumption.
   mr=DeformationModelReduction2DAxisymm



fens=FENodeSet(X);
fes=FESetQ4(conn=[1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
for ref=1:nref
    fens,fes=Q4refine(fens,fes);
    list=fenodeselect(fens,distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
    fens.xyz[list,:]= JFinEALE.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
end 
fes.axisymm=(mr== DeformationModelReduction2DAxisymm)  # note that this reflects the chosen model reduction

#     File =  "mesh.vtk"
# vtkexportmesh(File, fens, fes)

# now we create the geometry and displacement fields
geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

# Apply EBC's
l1=fenodeselect(fens,box=[-Inf Inf 0 0],inflate=tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])    
l1=fenodeselect(fens,box=[-Inf Inf 1.79  1.79],inflate=tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])    
applyebc!(u)
numberdofs!(u)

# Temperature field
dT =NodalField(name="dT",data=reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


# Property and material
material=MaterialDeformationLinear (PropertyDeformationLinearIso(E=Ea,nu=nua,CTE=alphaa))

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

K =stiffness(mr, femm, geom, u)
F = thermalstrainloads(mr, femm, geom, u, dT)
#K=cholfact(K)
U=  K\F
scattersysvec!(u,U[:])

nA =fenodeselect(fens,box=JFFlt[1.0  1.0 0.0 0.0],inflate=tolerance);

fld= fieldfromintegpoints(mr, femm, geom, u, dT, :Cauchy, 2)


# File =  "LE11NAFEMS_Q4_sigmay.vtk"
# vtkexportmesh(File, fens, fes; scalars=fld.values,scalars_name ="sigmay", vectors=u.values,vectors_name="u")

sA = fld.values[nA]/phun("MEGA*Pa")
sAn = fld.values[nA]/sigmaA
#println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
@test abs(sA[1]- (-86.1980522942295))<1.e-5

pt=1
function inspector(idat,out,loc,pc)
    global pt
    if (pt==1)
        @test norm(loc-[1.0419351909219474 0.010804307080807748])<1.e-7
        @test norm(reshape(out,1,4)-[2.293918786792204e6 -8.730762977944845e7 5.367202312538694e7 -1.1411272707766877e6])<1.e-2
        pt=pt+1
    end
    return idat
end

fen2fe =FENodeToFEMap(fes.conn,nfens(geom))
inspectintegpoints(mr,femm, geom, u, dT,  fen2fe.map[nA[1]],
                   inspector, []; output=:Cauchy)

#finealemesh(fens,fes,"meshmfile")

# end
# LE11NAFEMS()

end


