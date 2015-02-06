using JFinEALE
using JFinEALE.AlgoHeatDiffusionModule
using Base.Test

function Poisson1()
    t0 = time()

    A= 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity= eye(2,2); # conductivity matrix
    magn = reshape([-6.0],1,1); #heat source
    tempf(x)=(1.0 + x[:,1].^2 + 2*x[:,2].^2);#the exact distribution of temperature
    N=100;# number of subdivisions along the sides of the square domain.  Note that for the test problem N=1000.

    # println("""

    #         Heat conduction example described by Amuthan A. Ramabathiran
    #         http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    #         Unit square, with known temperature distribution along the boundary, 
    #         and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
    #         in a grid of $N x $N edges. 
    #         """
    #         )
    
    #println("Mesh generation")
    fens,fes =T3block(A, A, N, N)

    geom = NodalField(name ="geom",data =fens.xyz)
    Temp = NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))

    #println("Searching nodes  for BC")
    l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List=[l1, l2, l3, l4];
    setebc!(Temp,List,trues(length(List)),List*0+1,tempf(geom.values[List,:])[:])
    applyebc!(Temp)
    numberdofs!(Temp)

    t1 = time()

    p=PropertyHeatDiffusion(thermal_conductivity)
    material=MaterialHeatDiffusion (p)

    femm = FEMMHeatDiffusion(FEMMBase(fes, TriRule(npts=1)), material)


    #println("Conductivity")
    K=conductivity(femm, geom, Temp)
    #println("Nonzero EBC")
    F2=nzebcloadsconductivity(femm, geom, Temp);
    #println("Internal heat generation")
    fi = ForceIntensity(magn);
    F1=distribloads(femm.femmbase, geom, Temp, fi, 3);

    #println("Factorization")
    K=cholfact(K)
    #println("Solution of the factorized system")
    U=  K\(F1+F2)
    scattersysvec!(Temp,U[:])

    #println("Total time elapsed = $(time() - t0) [s]")
    #println("Solution time elapsed = $(time() - t1) [s]")

    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error=Error+abs(Temp.values[k,1]-tempf(fens.xyz[k,:]))
    end
    #println("Error =$Error")
    @test Error[1]<1.e-5

    true
end
Poisson1()

function annulus_Q4_example_algo()


    # #println("""

    #         Annular region, ingoing and outgoing flux. Temperature at one node prescribed. 
    #         Minimum/maximum temperature ~(+/-)0.591.
    #         Mesh of linear quadrilaterals.
    #         This version uses the JFinEALE algorithm module.
    #         """)

    t0 = time()

    kappa=0.2*[1.0 0.0; 0.0 1.0]; # conductivity matrix
    magn = 0.06;# heat flux along the boundary
    rin= 1.0;#internal radius
    rex= 2.0;#external radius
    nr=10; nc=80;
    Angle=2*pi;
    thickness= 1.0;
    tolerance=min(rin/nr, rin/nc/2/pi)/10000;



    fens,fes = Q4annulus(rin,rex,nr,nc,Angle)
    setotherdimension!(fes,thickness)
    fens,fes = mergenodes(fens, fes, tolerance);
    edge_fes = meshboundary (fes);
    setotherdimension!(edge_fes,thickness)

    # At a single point apply an essential boundary condition (pin down the temperature)
    l1 =fenodeselect(fens; box=[0.0 0.0 -rex -rex], inflate = tolerance)
    essential1=dmake(node_list=l1,temperature=0.0)

    # The flux boundary condition is applied at two pieces of surface
    # Side 1
    l1=feselect(fens,edge_fes,box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMBase(subset(edge_fes,l1),GaussRule(order=2,dim=1))
    flux1=dmake(femm=el1femm,normal_flux=-magn) # entering the domain
    # Side 2
    l2=feselect(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el2femm = FEMMBase(subset(edge_fes,l2),GaussRule(order=2,dim=1))
    flux2=dmake(femm=el2femm,normal_flux=+magn) # leaving the domain
    
    material=MaterialHeatDiffusion (PropertyHeatDiffusion(kappa))
    femm = FEMMHeatDiffusion(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

    # Make model data
    modeldata= dmake(fens= fens,
                     region=[dmake(femm=femm)],
                     boundary_conditions=dmake(flux=[flux1,flux2],essential=[essential1]));

    # Call the solver
    modeldata=AlgoHeatDiffusionModule.steadystate(modeldata)
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    #println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

    @test abs(minimum(Temp.values)-(-0.5911))<1.0e-4
    @test abs(maximum(Temp.values)-(+0.5911))<1.0e-4

    #println("Total time elapsed = ",time() - t0,"s")

    # Postprocessing
    # MeshExportModule.vtkexportmesh ("annulusmod.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")
end
annulus_Q4_example_algo()

function Poisson_FE_Q4_1()


    t0 = time()

    A= 1.0
    thermal_conductivity= eye(2,2); # conductivity matrix
    magn = -6.0; #heat source
    tempf(x)=(1.0 + x[:,1].^2 + 2*x[:,2].^2);
    N=100;

    # #println("""

    #         Heat conduction example described by Amuthan A. Ramabathiran
    #         http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    #         Unit square, with known temperature distribution along the boundary, 
    #         and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
    #         in a grid of $N x $N edges. 
    #         """
    #         )

    #println("Mesh generation")
    fens,fes =Q4block(A, A, N, N)



    geom = NodalField(name ="geom",data =fens.xyz)
    Temp = NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))


    #println("Searching nodes  for BC")
    l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List=[l1, l2, l3, l4];
    setebc!(Temp,List,trues(length(List)),List*0+1,tempf(geom.values[List,:])[:])
    applyebc!(Temp)

    numberdofs!(Temp)

    t1 = time()
    
    m=MaterialHeatDiffusion (PropertyHeatDiffusion(thermal_conductivity))
    femm = FEMMHeatDiffusion(FEMMBase(fes, GaussRule(order=2,dim=2)), m)

    #println("Conductivity")
    K=conductivity(femm, geom, Temp)
    #Profile.print()

    #println("Nonzero EBC")
    F2=nzebcloadsconductivity(femm, geom, Temp);
    #println("Internal heat generation")
    fi = ForceIntensity(magn);
    F1=distribloads(femm.femmbase, geom, Temp, fi, 3);


    #println("Factorization")
    K=cholfact(K)
    #println("Solution of the factorized system")
    U=  K\(F1+F2)
    scattersysvec!(Temp,U[:])


    #println("Total time elapsed = $(time() - t0) [s]")
    #println("Solution time elapsed = $(time() - t1) [s]")

    # using MeshExportModule

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error=Error+abs(Temp.values[k,1]-tempf(fens.xyz[k,:]))
    end
    #println("Error =$Error")
    @test Error[1]<1.e-5

    true
end
Poisson_FE_Q4_1()

function Poisson_FE_example_algo()


    A= 1.0
    thermal_conductivity=eye(2,2); # conductivity matrix
    magn = reshape([-6.0],1,1); #heat source
    boundaryf(x)=1.0 + x[:,1].^2 + 2.0*x[:,2].^2;
    N=20;

    # #println("""

    #         Heat conduction example described by Amuthan A. Ramabathiran
    #         http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    #         Unit square, with known temperature distribution along the boundary, 
    #         and uniform heat generation rate inside.  Mesh of regular TRIANGLES,
    #         in a grid of $N x $N edges. 
    #         This version uses the JFinEALE algorithm module.
    #         """
    #         )
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
                                   Q=(x,J,l)->magn,fes=fes,integration_rule=TriRule(npts=1))],
                     boundary_conditions=dmake(essential=[essential1]));


    # Call the solver
    modeldata=AlgoHeatDiffusionModule.steadystate(modeldata)

    #println("Total time elapsed = ",time() - t0,"s")

    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    femm=modeldata["region"][1]["femm"]
    function errfh(loc,val)
        exact = boundaryf(loc)
        return abs(((exact-val)*exact)[1])
    end

    femm.femmbase.integration_rule=TriRule(npts=6)
    E= integratefieldfunction (femm.femmbase,geom,Temp,errfh,0.0,m=3)
    #println("Error=$E")
    @test E<0.0025

    # Postprocessing
    # geom=modeldata["geom"]
    # Temp=modeldata["temp"]
    # MeshExportModule.vtkexportmesh ("a.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

end
Poisson_FE_example_algo()

function PoissonRm2()
    

    # println("""

    #         Heat conduction example described by Amuthan A. Ramabathiran
    #         http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    #         Unit square, with known temperature distribution along the boundary, 
    #         and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
    #         in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom). 

    #         With material orientation matrix supplied by a function.
    #             """
    #             )
    t0 = time()

    A= 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity= eye(2,2); # conductivity matrix
    magn = -6.0; #heat source
    tempf(x)=(1.0 + x[:,1].^2 + 2*x[:,2].^2);#the exact distribution of temperature
    N=100;# number of subdivisions along the sides of the square domain
    Rm=[-0.9917568452513019 -0.12813414805267656;    -0.12813414805267656 0.9917568452513019];
    #Rm=[-0.8020689950104449 -0.5972313850116512;    -0.5972313850116512 0.8020689950104447];
    function Rmfun(XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        return Rm
    end

    #println("Mesh generation")
    fens,fes =T3block(A, A, N, N)

    geom = NodalField(name ="geom",data =fens.xyz)
    Temp = NodalField(name ="Temp",data =zeros(size(fens.xyz,1),1))

    #println("Searching nodes  for BC")
    l1 =fenodeselect(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    l2 =fenodeselect(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    l3 =fenodeselect(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    l4 =fenodeselect(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List=[l1, l2, l3, l4];
    setebc!(Temp,List,trues(length(List)),List*0+1,tempf(geom.values[List,:])[:])
    applyebc!(Temp)
    numberdofs!(Temp)

    t1 = time()

    p=PropertyHeatDiffusion(thermal_conductivity)
    material=MaterialHeatDiffusion (p)

    femm = FEMMHeatDiffusion(FEMMBase(fes, TriRule(npts=1), MaterialOrientation(Rmfun)), material)


    #println("Conductivity")
    K=conductivity(femm, geom, Temp)
    #println("Nonzero EBC")
    F2=nzebcloadsconductivity(femm, geom, Temp);
    #println("Internal heat generation")
    fi = ForceIntensity(magn);
    F1=distribloads(femm.femmbase, geom, Temp, fi, 3);

    #println("Factorization")
    K=cholfact(K)
    #println("Solution of the factorized system")
    U=  K\(F1+F2)
    scattersysvec!(Temp,U[:])

    # println("Total time elapsed = $(time() - t0) [s]")
    # println("Solution time elapsed = $(time() - t1) [s]")

    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error=Error+abs(Temp.values[k,1]-tempf(fens.xyz[k,:]))
    end
    @test Error[1]<0.000025
    # println("Error =$Error")

    # using MeshExportModule

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

end
PoissonRm2();

module scratch
using JFinEALE
using Base.Test

## Two-dimensional heat transfer with convection: convergence study
#

## Description
#
# Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
# short edge the temperature is fixed at 100 °C, and on one long edge the
# plate is perfectly insulated so that the heat flux is zero through that
# edge. The other two edges are losing heat via convection to an ambient
# temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
# .°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
# There is no internal generation of heat. Calculate the temperature 0.2 m
# along the un-insulated long side, measured from the intersection with the
# fixed temperature side. The reference result is 18.25 °C.

##
# The reference temperature at the point A  is 18.25 °C according to the
# NAFEMS publication ( hich cites the book Carslaw, H.S. and J.C. Jaeger,
# Conduction of Heat in Solids. 1959: Oxford University Press).

##
# The present  tutorial will investigate the reference temperature  and it
# will attempt to  estimate the  limit value more precisely using a
# sequence of meshes and Richardson's extrapolation.

## Solution
#

kappa=[52 0; 0 52]*phun("W/(M*K)"); # conductivity matrix
h=750*phun("W/(M^2*K)");# surface heat transfer coefficient
Width=0.6*phun("M");# Geometrical dimensions
Height=1.0*phun("M");
HeightA=0.2*phun("M");
Thickness=1.0*phun("M");
tolerance =Width/1000;


# println("""
#         Two-dimensional heat transfer with convection: convergence study.
#         """
#         )
        t0 = time()

modeldata=[]
resultsTempA=JFFlt[]
for nref=1:5

 # The mesh is created from two quadrilaterals to begin with
fens,fes =Q4quadrilateral([0.0 0.0; Width 0.0; Width HeightA; 0.0 HeightA], 1,1)
fens2,fes2 =Q4quadrilateral([Width HeightA; Width Height; 0.0 Height; 0.0 HeightA], 1,1)
fens,newfes1,fes2= mergemeshes(fens, fes, fens2, fes2, tolerance)
fes=cat(newfes1,fes2)
# Refine the mesh desired number of times
for ref=1:nref
fens,fes=Q4refine(fens,fes);
end
fens,fes=Q4toT3(fens,fes);
bfes= meshboundary(fes)

    # Define boundary conditions

    ##
    # The prescribed temperature is applied along edge 1 (the bottom
    # edge in Figure 1)..

    l1 =fenodeselect(fens; box=[0. Width 0. 0.0], inflate =tolerance)
    essential1= dmake(node_list=l1,temperature=100);


    ##
    # The convection boundary condition is applied along the edges
    # 2,3,4. The elements along the boundary are quadratic line
    # elements L3. The order-four Gauss quadrature is sufficiently
    # accurate.
    l2 =feselect(fens,bfes; box=[Width Width  0.0 Height], inflate =tolerance)
    l3 =feselect(fens,bfes; box=[0.0 Width Height Height], inflate =tolerance)
    convection1= dmake(fes=subset(bfes,[l2,l3]),integration_rule=GaussRule(dim=1,order=2),
                       ambient_temperature=0,surface_transfer_coefficient=h);

# Make model data
modeldata= dmake(fens= fens,
region=[dmake(conductivity=kappa,
fes=fes,integration_rule=TriRule(npts=1))],
boundary_conditions=dmake(essential=[essential1],convection =[convection1]));

    # Call the solver
    modeldata=JFinEALE.AlgoHeatDiffusionModule.steadystate(modeldata)

    #println("Total time elapsed = ",time() - t0,"s")

    l4 =fenodeselect(fens; box=[Width Width HeightA HeightA], inflate =tolerance)


    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    

    ##
    # Collect the temperature  at the point A  [coordinates
    # (Width,HeightA)].
    push!(resultsTempA,Temp.values[l4][1]);

end

## 
# These are the computed results for the temperature at point A:
#println("$( resultsTempA  )")

@test (norm(resultsTempA-[22.787233619704548,19.181345788263776,18.51597493638965,18.38155620219601,18.306424697692535])
       )<1.0e-6
       
end

module Scratch

using JFinEALE
using Base.Test

## Two-dimensional heat transfer with convection: convergence study
#

## Description
#
# Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
# short edge the temperature is fixed at 100 °C, and on one long edge the
# plate is perfectly insulated so that the heat flux is zero through that
# edge. The other two edges are losing heat via convection to an ambient
# temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
# .°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
# There is no internal generation of heat. Calculate the temperature 0.2 m
# along the un-insulated long side, measured from the intersection with the
# fixed temperature side. The reference result is 18.25 °C.

##
# The reference temperature at the point A  is 18.25 °C according to the
# NAFEMS publication ( hich cites the book Carslaw, H.S. and J.C. Jaeger,
# Conduction of Heat in Solids. 1959: Oxford University Press).

##
# The present  tutorial will investigate the reference temperature  and it
# will attempt to  estimate the  limit value more precisely using a
# sequence of meshes and Richardson's extrapolation.

## Solution
#

kappa=[52 0; 0 52]*phun("W/(M*K)"); # conductivity matrix
h=750*phun("W/(M^2*K)");# surface heat transfer coefficient
Width=0.6*phun("M");# Geometrical dimensions
Height=1.0*phun("M");
HeightA=0.2*phun("M");
Thickness=1.0*phun("M");
tolerance =Width/1000;


# println("""
#         Two-dimensional heat transfer with convection: convergence study.
#         """
#         )
        t0 = time()

modeldata=[]
resultsTempA=JFFlt[]
for nref=1:3

    # The mesh is created from two triangles to begin with
    fens,fes =T3blockx ([0.0 ,Width] ,[0.0, HeightA])
    fens2,fes2 =T3blockx ([0.0 ,Width] ,[HeightA, Height])
    fens,newfes1,fes2= mergemeshes(fens, fes, fens2, fes2, tolerance)
    fes=cat(newfes1,fes2)
    # Refine the mesh desired number of times
    for ref=1:nref
        fens,fes=T3refine(fens,fes);
    end
    fens,fes=T3toT6(fens,fes);
    bfes= meshboundary(fes)

    # Define boundary conditions

    ##
    # The prescribed temperature is applied along edge 1 (the bottom
    # edge in Figure 1)..

    l1 =fenodeselect(fens; box=[0. Width 0. 0.0], inflate =tolerance)
    essential1= dmake(node_list=l1,temperature=100);


    ##
    # The convection boundary condition is applied along the edges
    # 2,3,4. The elements along the boundary are quadratic line
    # elements L3. The order-four Gauss quadrature is sufficiently
    # accurate.
    l2 =feselect(fens,bfes; box=[Width Width  0.0 Height], inflate =tolerance)
    l3 =feselect(fens,bfes; box=[0.0 Width Height Height], inflate =tolerance)
    convection1= dmake(fes=subset(bfes,[l2,l3]),integration_rule=GaussRule(dim=1,order=3),
                       ambient_temperature=0,surface_transfer_coefficient=h);

    # Make model data
    modeldata= dmake(fens= fens,
                     region=[dmake(conductivity=kappa,
                                   fes=fes,integration_rule=TriRule(npts=3))],
                     boundary_conditions=dmake(essential=[essential1],convection =[convection1]));


    # Call the solver
    modeldata=JFinEALE.AlgoHeatDiffusionModule.steadystate(modeldata)

    #println("Total time elapsed = ",time() - t0,"s")

    l4 =fenodeselect(fens; box=[Width Width HeightA HeightA], inflate =tolerance)


    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    

    ##
    # Collect the temperature  at the point A  [coordinates
    # (Width,HeightA)].
    push!(resultsTempA,Temp.values[l4][1]);

end

## 
# These are the computed results for the temperature at point A:
#println("$( resultsTempA  )")
@test (norm(resultsTempA- [17.902800196451892,18.33233979534389,18.29648316223162])       )<1.0e-6
      

end


