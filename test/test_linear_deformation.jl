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
    @test  abs(idat.s[1]- -910911.9777008593)<1.0e-9

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
