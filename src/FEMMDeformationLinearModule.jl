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


function nzebcloadsstiffness{MR<:DeformationModelReduction,
    S<:FESet}(::Type{MR},self::FEMMDeformationLinear{S}, 
              geom::NodalField{JFFlt},
              u::NodalField{JFFlt})
    assembler = SysvecAssembler()
    return  nzebcloadsstiffness(self,assembler,geom,u);
end
export nzebcloadsstiffness

function nzebcloadsstiffness{MR<:DeformationModelReduction,
    S<:FESet,A<:SysmatAssemblerBase}(::Type{MR},self::FEMMDeformationLinear{S}, assembler::A,
                                     geom::NodalField{JFFlt},
                                     u::NodalField{JFFlt})
    #  Compute load vector for nonzero EBC for fixed temperature.
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
    D::JFFltMat =zeros(JFFlt,nstr,nstr); # material stiffness matrix -- used as a buffer
    B::JFFltMat =zeros(JFFlt,nstr,Kedim); # strain-displacement matrix -- used as a buffer
    mo::MaterialOrientation=self.femmbase.mo;
    tangentmoduli!(MR,self.material,D;loc=[0.0]);# Moduli in material orientation
    startassembly!(assembler, Kedim, Kedim, nfes, u.nfreedofs, u.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(u,pu,conn);# retrieve element displacement vector
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
            gatherdofnumsasvec!(temp,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, -Ke*pu, dofnums); # assemble element load vector
        end
    end # Loop over elements
    F= makevector!(assembler);
    return F
end
export nzebcloadsstiffness

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
            # Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
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

function inspectintegpoints(self, 
        geom, u, dT, fe_list, context,
        inspector, idat)
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
   mo::MaterialOrientation=self.femmbase.mo;
     tangentmoduli!(MR,self.material,D;loc=[0.0]);# Moduli in material orientation
    # # Output orientation matrix
    # outputRm_identity = ~isfield(context,'outputRm');# if identity, work not needed
    # outputRm_constant = (outputRm_identity) ...
    #     || strcmp(class(context.outputRm),'double');# need to compute at each point
    # if (~outputRm_constant)
    #     outputRmh = context.outputRm;# handle to a function  to evaluate outputRm
    # else
    #     if (~outputRm_identity)
    #         outputRm = context.outputRm;
    #     end
    # end    
    # if ~isfield(context,'output')
    #     context.output=  'Cauchy';# default output
    # end
    # # Material
    # mat = self.material;
    # ms=[];
    # # Retrieve data for efficiency
    # conns = fes.conn; # connectivity
    # labels = fes.label; # finite element labels
    # xs =geom.values;
    # Us =u.values;
    # context.strain= [];
    # if isempty(dT)
    #     dTs=zeros(geom.nfens,1);
    # else
    #     dTs=dT.values;
    # end
    # # Now loop over selected fes 
    # for m=1:length(fe_list)
    #     i=fe_list(m);
    #     conn = conns(i,:); # connectivity
    #     conn =conns(i,:);
    #     x=xs(conn,:);
    #     U=reshape(u, gather_values(u,conn));
    #     dT =dTs(conn,:);
    #     # Loop over all integration points
    #     for j=1:npts
    #         c =Ns{j}'*x;# physical location
    #         u_c = transpose(Ns{j})*Us(conn,:);
    #         J = x' * Nders{j};
    #         if (~Rm_constant)# do I need to evaluate the local material orientation?
    #             if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
    #             else,                    Rm =Rmh(c,J,[]);                end
    #         end
    #         if (~outputRm_constant)# do I need to evaluate the output orientation?
    #             if (~isempty(labels )),  outputRm =outputRmh(c,J,labels(i));
    #             else,                    outputRm =outputRmh(c,J,[]);    end
    #         end
    #         Ndersp = Nders{j}/J;# derivatives wrt global coor
    #         Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
    #         B =self.hBlmat(self,Ns{j},Ndersp*Rm,c,Rm);# strains  in material coordinates, displacements in global coordinates
    #         context.strain = B*U;#strain in the material coords
    #         context.dT = transpose(Ns{j})*dT;
    #         context.xyz =c;
    #         context.ms=ms;
    #         [out,ignore] = update(mat, context.ms, context);
    #         switch context.output
    #             case 'Cauchy'
    #                 out =mat.stress_vector_rotation(Rm')*out;#  To global coordinate system
    #                 if (~outputRm_identity)
    #                     out =mat.stress_vector_rotation(outputRm)*out;# To output coordinate system
    #                 end
    #         end
    #         if ~isempty (inspector)
    #             idat =feval(inspector,idat,out,c,u_c,pc(j,:));
    #         end
    #     end
    # end
    return;
end


end
