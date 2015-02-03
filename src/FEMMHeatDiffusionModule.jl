module FEMMHeatDiffusionModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.MaterialOrientationModule
using JFinEALE.FEMMBaseModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.MaterialHeatDiffusionModule
using JFinEALE.AssemblyModule

# Class for heat diffusion finite element modeling machine.
type FEMMHeatDiffusion{S<:FESet}
    femmbase::FEMMBase{S} # base finite element modeling machine
    material::MaterialHeatDiffusion # material object
end
export FEMMHeatDiffusion

function conductivity{S<:FESet}(self::FEMMHeatDiffusion{S}, 
                                geom::NodalField{JFFlt}, temp::NodalField{JFFlt})
    assembler = SysmatAssemblerSparseSymm();
    return conductivity(self,assembler,geom,temp);
end
export conductivity

function conductivity{S<:FESet,A<:SysmatAssemblerBase}(self::FEMMHeatDiffusion{S}, assembler::A,
                                geom::NodalField{JFFlt}, temp::NodalField{JFFlt})
    #  Compute the conductivity matrix.
    #
    #  function K = conductivity(self, assembler, geom, temp)
    #
    #  Returns K as a matrix.
    #  Arguments
    #            self  = heat diffusion model
    #            assembler = descendent of the sysmat_assembler class
    #            geom=geometry field
    #            temp=temperature field
    # code_typed(FEMMHeatDiffusionModule.conductivity,(FEMMHeatDiffusionModule.FEMMHeatDiffusion,Any,NodalFieldModule.NodalField,Nodal
    # FieldModule.NodalField))

    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(temp); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Material
    # Note that the thermal conductivity matrix is in the
    # local  material orientation coordinates.
    kappa_bar =  self.material.property.thermal_conductivity;
    # Prepare assembler and temporaries
    Ke::JFFltMat =zeros(JFFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    kappa_bargradNT::JFFltMat=zeros(JFFlt,mdim,nne); # intermediate result -- used as a buffer
    startassembly!(assembler, Kedim, Kedim, nfes, temp.nfreedofs, temp.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(Ke, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            updateRm!(self.femmbase.mo,loc,J,getlabel(fes,i)); # Material orientation matrix
            At_mul_B!(RmTJ, self.femmbase.mo.Rm, J); # local Jacobian matrix 
            # gradient WRT material coordinates
            FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
            A_mul_Bt!(kappa_bargradNT,kappa_bar,gradN); # intermediate result
            @inbounds for nx=1:Kedim # Do: Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN' ; only the upper triangle
                @inbounds for mx=1:nx # only the upper triangle
                    @inbounds for px=1:sdim
                        Ke[mx,nx] = Ke[mx,nx] + gradN[mx,px]*(Jac*w[j])*kappa_bargradNT[px,nx]
                    end
                end
            end
            @inbounds for nx=1:Kedim #  Do: Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN' ;
                @inbounds for mx=nx+1:Kedim # complete the lower triangle
                    Ke[mx,nx] = Ke[nx,mx]
                end
            end
        end # Loop over quadrature points
        gatherdofnumsasvec!(temp,dofnums,conn);# retrieve degrees of freedom
        assemble!(assembler, Ke, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export conductivity

function nzebcloadsconductivity{S<:FESet}(self::FEMMHeatDiffusion{S}, 
                                          geom::NodalField{JFFlt},
                                          temp::NodalField{JFFlt})
    assembler = SysvecAssembler()
    return  nzebcloadsconductivity(self,assembler,geom,temp);
end
export nzebcloadsconductivity

function nzebcloadsconductivity{S<:FESet,A<:SysvecAssemblerBase}(self::FEMMHeatDiffusion{S}, assembler::A,
                                          geom::NodalField{JFFlt},
                                          temp::NodalField{JFFlt})
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
    const ndn::JFInt= ndofn(temp); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const Kedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Material
    # Note that the thermal conductivity matrix is in the
    # local  material orientation coordinates.
    kappa_bar =  self.material.property.thermal_conductivity;
    # Prepare assembler and temporaries
    Ke::JFFltMat =zeros(JFFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Kedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    RmTJ::JFFltMat =zeros(JFFlt,mdim,mdim); # intermediate result -- used as a buffer
    gradN::JFFltMat =zeros(JFFlt,nne,mdim); # intermediate result -- used as a buffer
    pT::JFFltMat=zeros(JFFlt,Kedim,1);
    startassembly!(assembler, temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasvec!(temp,pT,conn);# retrieve element coordinates
        if norm(pT) != 0     # Is the load nonzero?
            gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
            fill!(Ke, 0.0);
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
                updateRm!(self.femmbase.mo,loc,J,getlabel(fes,i)); # Material orientation matrix 
                At_mul_B!(RmTJ, self.femmbase.mo.Rm, J); # local Jacobian matrix 
                # gradient WRT material coordinates
                FESetModule.gradN!(fes,gradN,gradNparams[j],RmTJ);#Do: gradN = gradNparams[j]/RmTJ;
                @inbounds for nx=1:Kedim # Do: Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN' ;
                    @inbounds for kx=1:sdim
                        @inbounds for px=1:sdim
                            @inbounds for mx=1:Kedim
                                Ke[mx,nx] = Ke[mx,nx] + gradN[mx,px]*((Jac*w[j])*kappa_bar[px,kx])*gradN[nx,kx]
                            end
                        end
                    end
                end
            end # Loop over quadrature points
            gatherdofnumsasvec!(temp,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, -Ke*pT, dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end
export nzebcloadsconductivity

function surfacetransfer{S<:FESet}(self::FEMMHeatDiffusion{S}, 
                                   geom::NodalField{JFFlt},
                                   temp::NodalField{JFFlt},surfacetransfercoeff::JFFlt)
    assembler = SysmatAssemblerSparseSymm()
    return  surfacetransfer(self,assembler,geom,temp,surfacetransfercoeff);
end
export surfacetransfer

function surfacetransfer{S<:FESet,A<:SysmatAssemblerBase}(self::FEMMHeatDiffusion{S}, assembler::A,
                                                          geom::NodalField{JFFlt},
                                                          temp::NodalField{JFFlt},surfacetransfercoeff::JFFlt)
    # Compute the surface heat transfer matrix.
    #
    # function H = surface_transfer(self, assembler, geom, temp)
    #
    # Arguments
    #           self  = heat diffusion model  
    #           assembler = descendent of the sysmat_assembler class
    #           geom=geometry field
    #           temp=temperature field
    #
    # Returns H as a matrix.
    
    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(temp); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const Hedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    He::JFFltMat =zeros(JFFlt,Hedim,Hedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Hedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Hedim, Hedim, nfes, temp.nfreedofs, temp.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(He, 0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points 
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobiansurface(fes,conn, Ns[j], J, x);# Jacobian
            @inbounds for nx=1:Hedim # Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}' ; only the upper triangle
                @inbounds for mx=1:nx # only the upper triangle
                    He[mx,nx] = He[mx,nx] + (surfacetransfercoeff*Jac*w[j])*Ns[j][mx]*Ns[j][nx];
                end
            end
            @inbounds for nx=1:Hedim #  Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}'
                @inbounds for mx=nx+1:Hedim # complete the lower triangle
                    He[mx,nx] = He[nx,mx]
                end
            end
        end # Loop over quadrature points
        gatherdofnumsasvec!(temp,dofnums,conn);# retrieve degrees of freedom
        assemble!(assembler, He, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end
export surfacetransfer

function surfacetransferloads{S<:FESet}(self::FEMMHeatDiffusion{S}, 
                                        geom::NodalField{JFFlt},
                                        temp::NodalField{JFFlt},
                                        ambtemp::NodalField{JFFlt},
                                        surfacetransfercoeff::JFFlt)
    assembler = SysvecAssembler()
    return  surfacetransferloads(self,assembler,geom,temp,ambtemp,surfacetransfercoeff);
end
export surfacetransferloads

function surfacetransferloads{S<:FESet,
    A<:SysvecAssemblerBase}(self::FEMMHeatDiffusion{S}, assembler::A,
                            geom::NodalField{JFFlt},
                            temp::NodalField{JFFlt},
                            ambtemp::NodalField{JFFlt},
                            surfacetransfercoeff::JFFlt)
    # Compute the load vector corresponding to surface heat transfer.
#
# function F = surface_transfer_loads (self, assembler, geom, temp, amb)
#
#    Arguments
#     assembler =  descendent of sysvec_assembler
#     geom=geometry field
#     temp=temperature field
#     amb = ambient temperature field
#
 
    fes=self.femmbase.fes;                       # the finite element set
    # Constants
    const nfes::JFInt=count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(temp); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const Hedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    He::JFFltMat =zeros(JFFlt,Hedim,Hedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Hedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    pT::JFFltMat=zeros(JFFlt,Hedim,1);
    startassembly!(assembler, temp.nfreedofs);
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasvec!(ambtemp,pT,conn);# retrieve element coordinates
        if norm(pT) != 0     # Is the load nonzero?
            gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
            fill!(He, 0.0); # Initialize element matrix
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobiansurface(fes,conn, Ns[j], J, x);# Jacobian
                @inbounds for nx=1:Hedim # Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}' ; only the upper triangle
                    @inbounds for mx=1:nx # only the upper triangle
                        He[mx,nx] = He[mx,nx] + (surfacetransfercoeff*Jac*w[j])*Ns[j][mx]*Ns[j][nx];
                    end
                end
                @inbounds for nx=1:Hedim #  Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}'
                    @inbounds for mx=nx+1:Hedim # complete the lower triangle
                        He[mx,nx] = He[nx,mx]
                    end
                end
            end # Loop over quadrature points
            gatherdofnumsasvec!(temp,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, He*pT, dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end
export surfacetransferloads


function nzebcsurfacetransferloads{S<:FESet}(self::FEMMHeatDiffusion{S}, 
                                          geom::NodalField{JFFlt},
                                          temp::NodalField{JFFlt},
                            surfacetransfercoeff::JFFlt)
    assembler = SysvecAssembler()
    return  nzebcsurfacetransferloads(self,assembler,geom,temp,surfacetransfercoeff);
end
export nzebcsurfacetransferloads

function nzebcsurfacetransferloads{S<:FESet,A<:SysvecAssemblerBase}(self::FEMMHeatDiffusion{S}, assembler::A,
                                          geom::NodalField{JFFlt},
                                          temp::NodalField{JFFlt},
                            surfacetransfercoeff::JFFlt)
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
    const ndn::JFInt= ndofn(temp); # number of degrees of freedom per node
    const nne::JFInt =nfensperfe(fes); # number of nodes for element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes); # manifold dimension of the element
    const Hedim::JFInt =ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.femmbase);
    # Prepare assembler and temporaries
    He::JFFltMat =zeros(JFFlt,Hedim,Hedim);                # element matrix -- used as a buffer
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Hedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    pT::JFFltMat=zeros(JFFlt,Hedim,1);
    startassembly!(assembler, temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i=1:nfes # Loop over elements
        getconn!(fes,conn,i);
        gathervaluesasvec!(temp,pT,conn);# retrieve element coordinates
        if norm(pT) != 0     # Is the load nonzero?
            gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
            fill!(He, 0.0);
            for j=1:npts # Loop over quadrature points 
                At_mul_B!(loc,Ns[j],x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
                Jac = FESetModule.Jacobiansurface(fes,conn, Ns[j], J, x);# Jacobian
                @inbounds for nx=1:Hedim # Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}' ; only the upper triangle
                    @inbounds for mx=1:nx # only the upper triangle
                        He[mx,nx] = He[mx,nx] + (surfacetransfercoeff*Jac*w[j])*Ns[j][mx]*Ns[j][nx];
                    end
                end
                @inbounds for nx=1:Hedim #  Do: He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}'
                    @inbounds for mx=nx+1:Hedim # complete the lower triangle
                        He[mx,nx] = He[nx,mx]
                    end
                end
            end # Loop over quadrature points
            gatherdofnumsasvec!(temp,dofnums,conn); # retrieve degrees of freedom
            assemble!(assembler, He*pT, dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end
export nzebcsurfacetransferloads


end
