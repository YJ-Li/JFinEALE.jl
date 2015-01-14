module FEMMAcousticsModule

import Base.Complex

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FEMMBaseModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.MaterialAcousticFluidModule
using JFinEALE.AssemblyModule

# Class for linear acoustics finite element modeling machine.
type FEMMAcoustics{S<:FESet}
    femmbase::FEMMBase{S}
    material::MaterialAcousticFluid # material object
end
export FEMMAcoustics

function acousticmass{S<:FESet,T<:Number}(self::FEMMAcoustics{S}, 
                                          geom::NodalFieldModule.NodalField,
                                          P::NodalFieldModule.NodalField{T})
    # Make the default assembler object.
    assembler = SysmatAssemblerSparseSymm();
    return acousticmass(self,assembler,geom,P);
end

"""
Compute the acoustic mass matrix.

Return K as a matrix.
Arguments
self  = acoustics model
assembler = matrix assembler
geom=geometry field
P=acoustic (perturbation) pressure field
"""
function acousticmass{S<:FESet,T<:Number,A<:SysmatAssemblerBase}(self::FEMMAcoustics{S}, assembler::A,
                                                                 geom::NodalFieldModule.NodalField,
                                                                 P::NodalFieldModule.NodalField{T})

    fes=self.femmbase.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(P); # number of degrees of freedom per node
    const nne::JFInt = nfense(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    const Cedim::JFInt =ndn*nne;                          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    Ce::JFFltMat =zeros(JFFlt,Cedim,Cedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Cedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    Jinv::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    Jac::JFFlt=0.0;
    afactor::JFFlt=0.0;
    startassembly!(assembler, Cedim, Cedim, nfes, P.nfreedofs, P.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);# retrieve element node numbers
        gathervalues!(geom,x,conn);# retrieve element coordinates
        fill!(Ce, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            # gradient WRT global Cartesian coordinates
            FESetModule.gradN!(fes,gradN,gradNparams[j],J);#Do: gradN = gradNparams[j]/J;
            afactor=(Jac*w[j]);
            for nx=1:Cedim # Do: Ce = Ce + gradN*((Jac*w[j]))*gradN' ;
                for px=1:sdim
                    for mx=1:Cedim
                        Ce[mx,nx] = Ce[mx,nx] + gradN[mx,px]*(afactor)*gradN[nx,px]
                    end
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(P,dofnums,conn);# retrieve degrees of freedom
        assemble!(assembler, Ce, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export acousticmass

function nzebcloadsacousticmass{S<:FESet,T<:Number}(self::FEMMAcoustics{S},
                                                    geom::NodalField,
                                                    P::NodalField{T})
    assembler = SysvecAssembler(P.values[1])
    return nzebcloadsacousticmass(self,assembler,geom,P);
end
export nzebcloadsacousticmass

# """
# Compute load vector for nonzero EBC for fixed pressure.

#     %    Arguments
#     %      assembler =  descendent of sysvec_assembler
#     %      geom=geometry field
#     %      P=acoustic (perturbation) pressure field
#     %
#     % Return the assembled system vector F.
#     """
function nzebcloadsacousticmass{S<:FESet,T<:Number,A<:SysvecAssemblerBase}(self::FEMMAcoustics{S}, assembler::A,
                                                                           geom::NodalField,
                                                                           P::NodalField{T})
    fes=self.femmbase.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(P); # number of degrees of freedom per node
    const nne::JFInt = nfense(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= mdim(fes);     # manifold dimension of the element
    const Cedim::JFInt =ndn*nne;                          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = FEMMBaseModule.integrationdata(self.femmbase);
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    Ce::JFFltMat =zeros(JFFlt,Cedim,Cedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Cedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    pP::JFFltMat=zeros(JFFlt,Cedim,1);
    startassembly!(assembler, P.nfreedofs);
    # Now loop over all finite elements in the set
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);# retrieve element node numbers
        gathervalues!(P,pP,conn);# retrieve element coordinates
        if norm(pP) != 0     # Is the load nonzero?
            gathervalues!(geom,x,conn);# retrieve element coordinates
            fill!(Ce, 0.0);
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
                # gradient WRT global Cartesian coordinates
                FESetModule.gradN!(fes,gradN,gradNparams[j],J);#Do: gradN = gradNparams[j]/J;
                afactor::JFFlt=(Jac*w[j]);
                for nx=1:Cedim # Do: Ce = Ce + gradN*((Jac*w[j]))*gradN' ;
                    for px=1:sdim
                        for mx=1:Cedim
                            Ce[mx,nx] = Ce[mx,nx] + gradN[mx,px]*(afactor)*gradN[nx,px]
                        end
                    end
                end
            end # Loop over quadrature points
            gatherdofnums!(P,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, -Ce*pP, dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end
export nzebcloadsacousticmass


function acousticstiffness{S<:FESet,T<:Number}(self::FEMMAcoustics{S}, 
                                               geom::NodalFieldModule.NodalField,
                                               Pddot::NodalFieldModule.NodalField{T})
    # Make the default assembler object.
    assembler = SysmatAssemblerSparseSymm();
    return acousticstiffness(self,assembler,geom,Pddot);
end
export acousticstiffness


function acousticstiffness{S<:FESet,T<:Number,A<:SysmatAssemblerBase}(self::FEMMAcoustics{S}, assembler::A,
                                                                      geom::NodalFieldModule.NodalField,
                                                                      Pddot::NodalFieldModule.NodalField{T})
    #     """
    #     Compute the acoustic stiffness matrix.

# Return K as a matrix.
# Arguments
# self  = acoustics model
# assembler = matrix assembler
# geom=geometry field
# Pddot=acoustic (perturbation) second-derivative pressure field
# """
    
    fes=self.femmbase.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(Pddot); # number of degrees of freedom per node
    const nne::JFInt = nfense(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    const Sedim::JFInt = ndn*nne;               # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Material
    bulk_modulus =  self.material.property.bulk_modulus;
    mass_density =  self.material.property.mass_density;
    c = sqrt(bulk_modulus/mass_density); # sound speed
    oc2=1.0/c^2;
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    Se::JFFltMat =zeros(JFFlt,Sedim,Sedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Sedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Sedim, Sedim, nfes, Pddot.nfreedofs, Pddot.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);# retrieve element node numbers
        gathervalues!(geom,x,conn);# retrieve element coordinates
        fill!(Se, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            At_mul_B!(loc, Ns[j], x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            for nx=1:Sedim # Do: Se = Se + N*((Jac*w[j]))*N' ;
                for mx=1:Sedim
                    Se[mx,nx] = Se[mx,nx] + Jac*Ns[j][mx]*((oc2*w[j]))*Ns[j][nx]
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(Pddot,dofnums,conn);# retrieve degrees of freedom
        assemble!(assembler, Se, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export acousticstiffness

function nzebcloadsacousticstiffness{S<:FESet,T<:Number}(self::FEMMAcoustics{S},
                                                         geom::NodalFieldModule.NodalField,
                                                         Pddot::NodalFieldModule.NodalField{T})
    assembler = SysvecAssembler(Pddot.values[1])#T(0.0)
    return nzebcloadsacousticstiffness(self,assembler,geom,Pddot)
end
export nzebcloadsacousticstiffness

# """
# Compute load vector for nonzero EBC for fixed second-order pressure rate.

#     %    Arguments
#     %      assembler =  descendent of sysvec_assembler
#     %      geom=geometry field
#     %      Pddot=second order pressure rate
#     %
#     % Return the assembled system vector F.
#     """
function nzebcloadsacousticstiffness{S<:FESet,T<:Number,A<:SysvecAssemblerBase}(self::FEMMAcoustics{S}, assembler::A,
                                                                                geom::NodalFieldModule.NodalField,
                                                                                Pddot::NodalFieldModule.NodalField{T})
    fes=self.femmbase.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(Pddot); # number of degrees of freedom per node
    const nne::JFInt = nfense(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    const Sedim::JFInt =ndn*nne;                          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Material
    bulk_modulus =  self.material.property.bulk_modulus;
    mass_density =  self.material.property.mass_density;
    c = sqrt(bulk_modulus/mass_density); # sound speed
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfens);
    end
    # Prepare assembler and temporaries
    Se::JFFltMat =zeros(JFFlt,Sedim,Sedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Sedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    pPddot::JFMat{T}=zeros(T,Sedim,1);
    startassembly!(assembler, Pddot.nfreedofs);
    # Now loop over all finite elements in the set
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);# retrieve element node numbers
        gathervalues!(Pddot,pPddot,conn);# retrieve element coordinates
        if norm(pPddot) != 0     # Is the load nonzero?
            gathervalues!(geom,x,conn);# retrieve element coordinates
            fill!(Se, 0.0);
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
                for nx=1:Sedim # Do: Se = Se + gradN*((Jac*w[j]))*gradN' ;
                    for mx=1:Sedim
                        Se[mx,nx] = Se[mx,nx] + Ns[j][mx]*((Jac/c^2*w[j]))*Ns[j][nx]
                    end
                end
            end # Loop over quadrature points
            gatherdofnums!(Pddot,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, -Se*pPddot, dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end
export nzebcloadsacousticstiffness

function acousticABC{S<:FESet,T<:Number}(self::FEMMAcoustics{S}, 
                                         geom::NodalFieldModule.NodalField,
                                         Pdot::NodalFieldModule.NodalField{T})
    # Were we supplied assembler object?  If not make a default.
    assembler = SysmatAssemblerSparseSymm();
    return acousticABC(self,assembler,geom,Pdot);
end
export acousticABC

# """
#     Compute the acoustic ABC (Absorbing Boundary Condition) matrix.

#     Arguments
#     geom=geometry field
#     Pdot=time derivative of the acoustic pressure field
#     """
function acousticABC{S<:FESet,T<:Number,A<:SysmatAssemblerBase}(self::FEMMAcoustics{S}, assembler::A,
                                                                geom::NodalFieldModule.NodalField,
                                                                Pdot::NodalFieldModule.NodalField{T})


    fes=self.femmbase.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(Pdot); # number of degrees of freedom per node
    const nne::JFInt = nfense(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    const Dedim::JFInt =ndn*nne;                          # dimension of the element matrix
    # Were we supplied assembler object?  If not make a default.
    if assembler== nothing
        assembler = SysmatAssemblerSparseSymm()
    end
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Material
    bulk_modulus =  self.material.property.bulk_modulus;
    mass_density =  self.material.property.mass_density;
    c = sqrt(bulk_modulus/mass_density); # sound speed
    # Prepare some data:
    labels=fes.label; # individual element labels
    if length(fes.label)==0
        labels=zeros(JFInt,nfes);
    end
    # Prepare assembler and temporaries
    De::JFFltMat =zeros(JFFlt,Dedim,Dedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Dedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Dedim, Dedim, nfes, Pdot.nfreedofs, Pdot.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);# retrieve element node numbers
        gathervalues!(geom,x,conn);# retrieve element coordinates
        fill!(De, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobiansurface(fes,conn, Ns[j], J, x);# Jacobian
            for nx=1:Dedim # Do: De = De + N*((Jac*w[j]))*N' ;
                for mx=1:Dedim
                    De[mx,nx] = De[mx,nx] + Ns[j][mx]*((Jac/c*w[j]))*Ns[j][nx]
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(Pdot,dofnums,conn);# retrieve degrees of freedom
        assemble!(assembler, De, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export acousticABC

end
