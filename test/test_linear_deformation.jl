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


