module FEMMDeformationLinearModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FEMMBaseModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.AssemblyModule
using JFinEALE.DeformationModelReductionModule
using JFinEALE.PropertyDeformationLinearModule
using JFinEALE.MaterialDeformationLinearModule

# Class for heat diffusion finite element modeling machine.
type FEMMDeformationLinear{S<:FESet,P<:PropertyDeformationLinear}
    femmbase::FEMMBase{S} # base finite element modeling machine
    #modelreduction::DeformationModelReduction # deformation model reduction type
    material::MaterialDeformationLinear{P} # material object
end
export FEMMDeformationLinear

function stiffness{MR<:DeformationModelReduction,
    S<:FESet}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{JFFlt})
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(MR,self,assembler,geom,u);
end
export stiffness

function stiffness{MR<:DeformationModelReduction,
    S<:FESet,A<:SysmatAssemblerBase}(::Type{MR},
                                     self::FEMMDeformationLinear{S},
                                     assembler::A,
                                     geom::NodalField{JFFlt}, u::NodalField{JFFlt})
    #  Compute the stiffness matrix.
    #

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfense(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    Ke::JFFltMat =zeros(JFFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    Rm::JFFltMat=eye(JFFlt,sdim,mdim); # material orientation matrix -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    D::JFFltMat =zeros(JFFlt,nstr,nstr); # material stiffness matrix -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
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
            updateRm!(self.femmbase,loc,J,labels[i]); # Material orientation matrix 
            At_mul_B!(RmTJ, self.femmbase.Rm, J); # local Jacobian matrix 
            # gradient WRT material coordinates
            FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
            Blmat!(MR,B,Ns[j],gradN,loc,Rm);#  strains in material cs, displacements in global cs
            #tangentmoduli!(MR,mat,D,loc);# Moduli in material orientation
            for nx=1:Kedim # Do: Ke = Ke + (B'*(D*(Jac*w[j]))*B); 
                for kx=1:nstr
                    for px=1:nstr
                        for mx=1:Kedim
                            Ke[mx,nx] = Ke[mx,nx] + B[px,mx]*((Jac*w[j])*D[px,kx])*B[kx,nx]
                        end
                    end
                end
            end
        end # Loop over quadrature points
        gatherdofnumsasvec!(u, dofnums, conn);# retrieve degrees of freedom
        assemble!(assembler, Ke, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export stiffness


function mass{MR<:DeformationModelReduction,
    S<:FESet}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt}, u::NodalField{JFFlt})
    assembler = SysmatAssemblerSparseSymm();
    return mass(MR,self,assembler,geom,u);
end
export mass

function mass{MR<:DeformationModelReduction,
    S<:FESet,A<:SysmatAssemblerBase}(::Type{MR},
                                     self::FEMMDeformationLinear{S},
                                     assembler::A,
                                     geom::NodalField{JFFlt}, u::NodalField{JFFlt})
    #  Compute the mass matrix.
    #

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfense(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const nstr::JFInt =nstrains(MR);  # number of strains
    const Medim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    Me::JFFltMat =zeros(JFFlt,Medim,Medim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Medim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    Rm::JFFltMat=eye(JFFlt,sdim,mdim); # material orientation matrix -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    rho::JFFltMat =self.material.property.mass_density; # mass density
    NexpTNexp=Array(JFFltMat,1,npts);# basis f. matrix -- used as a buffer
    for j=1:npts # This quantity is the same for all quadrature points
        Nexp= zeros(JFFlt,ndn,Medim)
        for l1 = 1:nne
            Nexp[1:ndn,(l1-1)*ndn+1:(l1)*ndn]=eye(dim)*Ns{j}[l];
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
            Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
            thefactor::JFFlt =(rho*Jac*w[j]);
            for nx=1:Medim # Do: Me = Me + (B'*(D*(Jac*w[j]))*B); 
                for mx=1:Medim
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

end
