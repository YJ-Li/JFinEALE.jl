module FEMMDeformationWinklerModule

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
type FEMMDeformationWinkler{S<:FESet}
    femmbase::FEMMBase{S} # base finite element modeling machine
end
export FEMMDeformationWinkler



function surfacenormalspringstiffness{S<:FESet,T<:Number}(self::FEMMDeformationWinkler{S}, 
              geom::NodalField{JFFlt}, u::NodalField{T},springconstant::JFFlt)
    assembler = SysmatAssemblerSparseSymm();
    return surfacenormalspringstiffness(self,assembler,geom,u,springconstant);
end
export surfacenormalspringstiffness

function  surfacenormalspringstiffness{S<:FESet,
    A<:SysmatAssemblerBase,T<:Number}(self::FEMMDeformationWinkler{S},assembler::A, 
              geom::NodalField{JFFlt}, u::NodalField{T},springconstant::JFFlt)
# Compute the stiffness matrix of surface normal spring.
#
# function H = surface_normal_spring_stiffness(self, assembler, geom, u)
#
# Arguments
#           self  = heat diffusion model  
#           assembler = descendent of the sysmat_assembler class
#           geom=geometry field
#           u=displacement field
#
# Returns H as a matrix.
#
# Rationale: consider continuously distributed springs between the surface of 
# the solid body and the 'ground', in the direction normal to the surface. 
# If the spring coefficient becomes large, we have an approximate
    # method of enforcing the normal displacement to the surface.
    
    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(u); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
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
    startassembly!(assembler, Kedim, Kedim, nfes, u.nfreedofs, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(Ke, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobiansurface(fes,conn, Ns[j], J, x);# Jacobian
            n= surfacenormal(loc,J);# find the normal to the surface
            Nn =reshape(n*Ns[j]',Kedim,1);# The normal n is a column vector
            @inbounds for nx=1:Kedim # Do: He = He + Nn*(Nn'*(springconstant*(Jac*w(j))); only the upper triangle
                @inbounds for mx=1:nx # only the upper triangle
                    Ke[mx,nx] = Ke[mx,nx] + Nn[mx]*(springconstant*Jac*w[j])*Nn[nx]
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
export surfacenormalspringstiffness

function  surfacenormal(loc::JFFltMat,J::JFFltMat)
    # Compute local normal. This makes sense for bounding surfaces only.
    #
    # function norml = normal(self,c,J)
    #
    #      self=finite element block
    #      c=spatial location
    #      J= Jacobian matrix
    #

    norml= zeros(JFFlt, length(loc))
    # Produce a default normal
    if (size(J,1)==3)&&(size(J,2)==2)# surface in three dimensions
        norml= cross(vec(J[:,1]),vec(J[:,2]));# outer normal to the surface
        norml=norml/norm(norml);
    elseif (size(J,1)==2)  && (size(J,2)==1)# curve in two dimensions
        norml= [J[2,1];-J[1,1]];# outer normal to the contour
        norml=norml/norm(norml);
    else
        error("No definition of normal vector");
    end
    
    return reshape(norml,length(norml),1) # return a column vector
end


end
