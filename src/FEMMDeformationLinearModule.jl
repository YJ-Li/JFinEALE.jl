module FEMMDeformationLinearModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.MaterialOrientationModule
using JFinEALE.FEMMBaseModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.AssemblyModule
using JFinEALE.DeformationModelReductionModule
using JFinEALE.PropertyDeformationLinearModule
using JFinEALE.MaterialDeformationModule
using JFinEALE.MaterialDeformationLinearModule
using JFinEALE.FENodeToFEMapModule

# Class for heat diffusion finite element modeling machine.
type FEMMDeformationLinear{S<:FESet,P<:PropertyDeformationLinear}
    femmbase::FEMMBase{S} # base finite element modeling machine
    #modelreduction::DeformationModelReduction # deformation model reduction type
    material::MaterialDeformationLinear{P} # material object
end
export FEMMDeformationLinear

function qpadd!(Ke::JFFltMat,Kedim::JFInt,nstr::JFInt,B::JFFltMat,Jacw::JFFlt,D::JFFltMat)
    @inbounds for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B); 
        @inbounds for kx=1:nstr
            @inbounds for px=1:nstr
                @inbounds for mx=1:Kedim
                    Ke[mx,nx] = Ke[mx,nx] + B[px,mx]*((Jacw)*D[px,kx])*B[kx,nx]
                end
            end
        end
    end
end

function f()
    Ke=zeros(JFFlt,60,60)
    B=rand(JFFlt,6,60)
    Jac= 1.3375::JFFlt
    w= 0.211::JFFlt
    D=rand(JFFlt,6,6)
    Kedim=size(Ke,1)
    nstr=6
    @inbounds for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B); 
        @inbounds for kx=1:nstr
            @inbounds for px=1:nstr
                @inbounds for mx=1:Kedim
                    Ke[mx,nx] = Ke[mx,nx] + B[px,mx]*((Jac*w)*D[px,kx])*B[kx,nx]
                end
            end
        end
    end
end

function stiffness{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T})
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(MR,self,assembler,geom,u);
end
export stiffness

# @code_warntype stiffness(DeformationModelReduction,femm,SysmatAssemblerSparseSymm(),geom, u)
function stiffness{MR<:DeformationModelReduction,
    S<:FESet,A<:SysmatAssemblerBase,T<:Number}(::Type{MR},
                                     self::FEMMDeformationLinear{S},
                                     assembler::A,
                                     geom::NodalField{JFFlt}, u::NodalField{T})
    #  Compute the stiffness matrix.
    #

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    Ke::JFFltMat =zeros(JFFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    D::JFFltMat =zeros(JFFlt,nstr,nstr); # material stiffness matrix -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    DB::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    mo::MaterialOrientation=self.femmbase.mo;
    tangentmoduli!(MR,self.material,D;loc=[0.0]);# Moduli in material orientation
    startassembly!(assembler, Kedim, Kedim, nfes, u.nfreedofs, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(Ke, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            updateRm!(mo,loc,J,getlabel(fes,i)); # Material orientation matrix 
            At_mul_B!(RmTJ, mo.Rm, J); # local Jacobian matrix 
            # gradient WRT material coordinates
            FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
            Blmat!(MR,B,Ns[j],gradN,loc,mo.Rm);#  strains in material cs, displacements in global cs
            A_mul_B!(DB,D,B)    # intermediate product
            #tangentmoduli!(MR,mat,D,loc);# Moduli in material orientation
            @inbounds for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B); only the upper triangle
                @inbounds for mx=1:nx # only the upper triangle
                    @inbounds for px=1:nstr
                        Ke[mx,nx] = Ke[mx,nx] + B[px,mx]*(Jac*w[j])*DB[px,nx]
                    end
                end
            end
            @inbounds for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B);
                @inbounds for mx=nx+1:Kedim # complete the lower triangle
                    Ke[mx,nx] = Ke[nx,mx]
                end
            end 
        end # Loop over quadrature points
        gatherdofnumsasvec!(u, dofnums, conn);# retrieve degrees of freedom
        assemble!(assembler, Ke, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export stiffness


function nzebcloadsstiffness{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt},
              u::NodalField{T})
    assembler = SysvecAssembler()
    return  nzebcloadsstiffness(MR,self,assembler,geom,u);
end
export nzebcloadsstiffness

function nzebcloadsstiffness{MR<:DeformationModelReduction,
    S<:FESet,A<:SysvecAssemblerBase,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, assembler::A,
                                     geom::NodalField{JFFlt},
                                     u::NodalField{T})
    #  Compute load vector for nonzero EBC for fixed displacement.
    # 
    # %    Arguments
    # %      assembler =  descendent of sysvec_assembler
    # %      geom=geometry field
    # %      temp=temperature field
    # %
    # % Return the assembled system vector F.

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    Ke::JFFltMat =zeros(JFFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    pu::JFFltMat=zeros(JFFlt,Kedim,1);
    D::JFFltMat =zeros(JFFlt,nstr,nstr); # material stiffness matrix -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    mo::MaterialOrientation=self.femmbase.mo;
    tangentmoduli!(MR,self.material,D;loc=[0.0]);# Moduli in material orientation
    startassembly!(assembler, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasvec!(u,pu,conn);# retrieve element displacement vector
        if norm(pu) != 0     # Is the load nonzero?
            gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
            fill!(Ke, 0.0); # Initialize element matrix
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
                updateRm!(mo,loc,J,getlabel(fes,i)); # Material orientation matrix 
                At_mul_B!(RmTJ, mo.Rm, J); # local Jacobian matrix 
                # gradient WRT material coordinates
                FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
                Blmat!(MR,B,Ns[j],gradN,loc,mo.Rm);#  strains in material cs, displacements in global cs
                #tangentmoduli!(MR,mat,D,loc);# Moduli in material orientation
                @inbounds for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B); 
                    @inbounds for kx=1:nstr
                        @inbounds for px=1:nstr
                            @inbounds for mx=1:Kedim
                                Ke[mx,nx] = Ke[mx,nx] + B[px,mx]*((Jac*w[j])*D[px,kx])*B[kx,nx]
                            end
                        end
                    end
                end
            end # Loop over quadrature points
            gatherdofnumsasvec!(u,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, -Ke*pu, dofnums); # assemble element load vector
        end
    end # Loop over elements
    F= makevector!(assembler);
    return F
end
export nzebcloadsstiffness

function mass{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T})
    assembler = SysmatAssemblerSparseSymm();
    return mass(MR,self,assembler,geom,u);
end
export mass

function mass{MR<:DeformationModelReduction,
    S<:FESet,A<:SysmatAssemblerBase,T<:Number}(::Type{MR},
                                     self::FEMMDeformationLinear{S},
                                     assembler::A,
                                     geom::NodalField{JFFlt}, u::NodalField{T})
    #  Compute the mass matrix.
    #

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Medim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    Me::JFFltMat =zeros(JFFlt,Medim,Medim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Medim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    rho::JFFlt =self.material.property.mass_density; # mass density
    NexpTNexp=Array(JFFltMat,1,npts);# basis f. matrix -- used as a buffer
    for j=1:npts # This quantity is the same for all quadrature points
        Nexp= zeros(JFFlt,ndn,Medim)
        for l1 = 1:nne
            Nexp[1:ndn,(l1-1)*ndn+1:(l1)*ndn]=eye(ndn)*Ns[j][l1];
        end
        NexpTNexp[j]=Nexp'*Nexp;
    end
    startassembly!(assembler, Medim, Medim, nfes, u.nfreedofs, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(Me, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            thefactor::JFFlt =(rho*Jac*w[j]);
            @inbounds for nx=1:Medim # Do: Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
                @inbounds for mx=1:Medim
                    Me[mx,nx] = Me[mx,nx] + NexpTNexp[j][mx,nx]*thefactor
                end
            end
        end # Loop over quadrature points
        gatherdofnumsasvec!(u, dofnums, conn);# retrieve degrees of freedom
        assemble!(assembler, Me, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export mass


function thermalstrainloads{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt},
              u::NodalField{T},dT::NodalField{JFFlt})
    assembler = SysvecAssembler()
    return  thermalstrainloads(MR,self,assembler,geom,u,dT);
end
export thermalstrainloads

function  thermalstrainloads{MR<:DeformationModelReduction,
    S<:FESet,A<:SysvecAssemblerBase,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S},
                                               assembler::A,
                                               geom::NodalField{JFFlt},
                                               u::NodalField{T},dT::NodalField{JFFlt})
# Compute the thermal-strain load vectors of the finite element set.
#
# function F = thermal_strain_loads(self, assembler, geom, u, dT)
#
# Return the assembled system vector F.
#    Arguments
#     assembler =  descendent of sysvec_assembler
#     geom=geometry field
#     u=displacement field
#     dT=temperature difference field (current temperature minus the
#         reference temperature at which the solid experiences no strains)
  
    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    Fe::JFFltMat =zeros(JFFlt,Kedim,1);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    DeltaT::JFFltMat=zeros(JFFlt,nne,1) # array of nodal temperatures -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    pu::JFFltMat=zeros(JFFlt,Kedim,1);
    D::JFFltMat =zeros(JFFlt,nstr,nstr); # material stiffness matrix -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    mo::MaterialOrientation=self.femmbase.mo;
    tangentmoduli!(MR,self.material,D;loc=[0.0]);# Moduli in material orientation
    startassembly!(assembler, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasvec!(dT,DeltaT,conn);# retrieve element temperatures
        if norm(DeltaT) != 0     # Is the thermal increment nonzero?
            gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
            fill!(Fe, 0.0); # Initialize element matrix
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
                updateRm!(mo,loc,J,getlabel(fes,i)); # Material orientation matrix 
                At_mul_B!(RmTJ, mo.Rm, J); # local Jacobian matrix 
                # gradient WRT material coordinates
                FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
                Blmat!(MR,B,Ns[j],gradN,loc,mo.Rm);#  strains in material cs, displacements in global cs
                #tangentmoduli!(MR,mat,D,loc);# Moduli in material orientation
                iSigma=thermalstress(MR,self.material,dT= Ns[j]'*DeltaT)
                @inbounds for mx=1:Kedim
                    @inbounds for px=1:nstr
                        Fe[mx] = Fe[mx] - B[px,mx]*(Jac*w[j])*iSigma[px]
                    end
                end
            end
            gatherdofnumsasvec!(u,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, Fe, dofnums); # assemble element load vector
        end
    end # Loop over elements
    F= makevector!(assembler);
    return F
end
export thermalstrainloads


function inspectintegpoints{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T},
              felist::JFIntVec,
              inspector::Function, idat; context...)
    dT=NodalField(name ="dT",data =zeros(JFFlt,nfens(geom),1)) # zero difference in temperature
    return inspectintegpoints(MR,self,geom,u,dT,felist,inspector,idat; context...);
end
export inspectintegpoints

function inspectintegpoints{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T},
              dT::NodalField{JFFlt},
              felist::JFIntVec,
              inspector::Function, idat; context...)
    # Inspect the integration point quantities.
# 
# Input arguments
#    geom - reference geometry field
#    u - displacement field
#    dT - temperature difference field
#    fe_list - indexes of the finite elements that are to be inspected:
#          The fes to be included are: fes(fe_list).
#    context    - structure: see the update() method of the material.
#    inspector - function handle or in-line function with the signature
#             idat =inspector(idat, out, xyz, u, pc),
#        where
#         idat - a structure or an array that the inspector may
#                use to maintain some state, for instance minimum or
#                maximum of stress, out is the output  of the update()
#                method, xyz is the location of the integration point
#                in the *reference* configuration, and u the displacement 
#                of the integration point. The argument pc are the
#                parametric coordinates of the quadrature point.
# Output arguments
#     idat - see the description of the input argument

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
     conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    ue::JFFltMat =zeros(JFFlt,Kedim,1); # array of node displacements -- used as a buffer
    dTe::JFFltMat =zeros(JFFlt,nne,1); # array of node temperature increments -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    ms =[]                                # stand-in for a GENUINE material state
    mo::MaterialOrientation=self.femmbase.mo; # the material coordinate system
    outputmo=mo;                # output the stresses in the material coordinate system
    output=:Cauchy;
    for arg in context
        sy, val = arg
        if sy==:outputmo
            outputmo=val
        elseif sy==:output
            output=val
        end
    end
    
    for ilist=1:length(felist) # Loop over elements
        i=felist[ilist];
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        gathervaluesasvec!(u,ue,conn);# retrieve element displacements
        gathervaluesasvec!(dT,dTe,conn);# retrieve element temperature increments
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            updateRm!(mo,loc,J,getlabel(fes,i)); # Material orientation matrix 
            At_mul_B!(RmTJ, mo.Rm, J); # local Jacobian matrix
             # gradient WRT material coordinates
            FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
            Blmat!(MR,B,Ns[j],gradN,loc,mo.Rm);#  strains in material cs, displacements in global cs
            updateRm!(outputmo,loc,J,getlabel(fes,i)); # Output material orientation matrix
            # Quadrature point quantities
            qpstrain=B*ue;
            qpdT=Ns[j]'*dTe;
            # Material updates the state and returns the output
            out,newms =update!(MR,self.material,ms;
                               strain=qpstrain,dT=qpdT,loc=loc);
            if (output==:Cauchy)   # Transform stress tensor, if that is the output
                out =stressvectorrotation(MR,mo.Rm')*out;#  To global coordinate system
                out =stressvectorrotation(MR,outputmo.Rm)*out;# To output coordinate system
            end
            # Call the inspector
            idat =inspector(idat,out,loc,pc[j,:]); 
        end # Loop over quadrature points
    end # Loop over elements
    return idat;
end
export inspectintegpoints

function fieldfromintegpoints{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T},
              output, component; context...)
    dT=NodalField(name ="dT",data =zeros(JFFlt,nfens(geom),1)) # zero difference in temperature
    return fieldfromintegpoints(MR,self,geom,u,dT, output, component; context...)
end
export fieldfromintegpoints

function fieldfromintegpoints{MR<:DeformationModelReduction,
    S<:FESet,T<:Number}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T},
              dT::NodalField{JFFlt}, output, component; context...)
    
# Input arguments
#     geom     - reference geometry field
#     u        - displacement field
#     dT       - temperature difference field
#     output   - this is what you would assign to the 'output' attribute of
#                the context argument of the material update() method.
#     component- component of the 'output' array: see the material update()
#                method.
#     options  - struct with fields recognized by update() method of the 
#                material; optional argument
# Output argument
#     fld - the new field that can be used to map values the colors
#
   fes=self.femmbase.fes;
    # Container of intermediate results
    sum_inv_dist =zeros(JFFlt,nfens(geom));
    sum_quant_inv_dist =zeros(JFFlt,nfens(geom),length(component));
    # The data for the field to be constructed is initialized
    nvals = zeros(JFFlt,nfens(geom),length(component));
     # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Prepare temporaries
     conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    # This is an inverse-distance interpolation inspector.
    function idi(idat, out, xyz, pc)
        d=zeros(JFFlt,length(conn))
        mindn=Inf
        for jjj=1:length(d)
            d[jjj] =sum((x[jjj,:]-xyz).^2);
            if (d[jjj] > 0.0)
                mindn=min(mindn,d[jjj])
            end
        end
        mindn=mindn/1.0e9;
        for jjj=1:length(d)
            invdjjj =1./(d[jjj]+mindn);
            quant= out[component]
            for kkk=1:length(quant)
                sum_quant_inv_dist[conn[jjj],kkk]+= invdjjj*quant[kkk];
            end
            sum_inv_dist[conn[jjj]]=sum_inv_dist[conn[jjj]]+invdjjj;
        end
        return idat
    end
    # Loop over cells to interpolate to nodes
    idat=0;
    for i=1:count(fes)
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        idat = inspectintegpoints(MR, self, geom, u, dT,
                                  [i], idi, idat; output = output, context...);
    end
    # compute the data array
    for kkk=1:size(nvals,1)
        for j=1:length(component)
            if (sum_inv_dist[kkk]> 0.0)
                nvals[kkk,j]=sum_quant_inv_dist[kkk,j]/sum_inv_dist[kkk];
            end
        end
    end    
    # Make the field
    fld =NodalField(data=nvals);
    return fld;
end
export fieldfromintegpoints



end
