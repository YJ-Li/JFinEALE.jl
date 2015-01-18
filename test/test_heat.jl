using JFinEALE
using JFinEALE.HeatDiffusionAlgorithmModule
using Base.Test

function Poisson1()
    t0 = time()

    A= 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity= eye(2,2); # conductivity matrix
    magn = -6.0; #heat source
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
    modeldata=HeatDiffusionAlgorithmModule.steadystate(modeldata)
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
    magn = -6.0; #heat source
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
                                   Q=(x,J)->[magn],fes=fes,integration_rule=TriRule(npts=1))],
                     boundary_conditions=dmake(essential=[essential1]));


    # Call the solver
    modeldata=HeatDiffusionAlgorithmModule.steadystate(modeldata)

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
    function Rmfun!(XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
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

    femm = FEMMHeatDiffusion(FEMMBase(fes, TriRule(npts=1), Rmfun!), material)


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
