module FEMMBaseModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.IntegRuleModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.AssemblyModule
using JFinEALE.MaterialOrientationModule

#     % Base class for all finite element models
#     % from which more specialized (heat diffusion, elasticity,...) finite
#     % element model classes are derived.
#  Parameters:
#   mandatory:
#       fes = finite element set.  The type of the FE set will be dependent upon
#         the operations required. For instance, for interior (volume) integrals
#         such as body load or the stiffness hexahedral H8 may be used whereas
#         for boundary  (surface) integrals quadrilateral Q4 would be needed.
#   optional:
#       integration_rule= integration rule object
#       Rm= material orientation matrix, or a handle to function to compute the
#          material orientation matrix, or a string denoting the type
#          of material orientation matrix to be used ('geniso').
#             In the columns of the material orientation matrix are the basis vectors
#          expressed in the global Cartesian basis. If the orientation matrix
#          is not supplied, Rm==identity is assumed. The orientation matrix
#          may be passed in as empty ([]), and the model code can take advantage
#          of that by assuming that Rm being empty means the local-to-global
#          transformation is the identity and avoiding thereby the
#          multiplication.
#             The function to compute the orientation matrix should have
#          the signature
#                function Rm = SampleRm(XYZ, ts, label)
#          The orientation matrix can be computed based on any of the three
#          arguments.
#          XYZ= global Cartesian location of the point at which the orientation
#                  matrix is desired,
#          ts= the Jacobian matrix with the tangents to the parametric coordinates
#                  as columns,
    #          label= the label of the finite element in which the point XYZ resides

type FEMMBase{T<:FESet}
    fes::T # finite element set object
    integration_rule::IntegRule  # integration rule object
    mo::MaterialOrientation # updater of the material orientation matrix
    fen2fe_map  # map from finite element nodes to connected finite elements
    
    function  FEMMBase (fes,
                        integration_rule,
                        mo)
        return  new (deepcopy(fes), deepcopy(integration_rule), mo)
    end

end
export FEMMBase

# Construct with the user-defined orientation matrix updater.
function  FEMMBase{T<:FESet} (fes::T,
                              integration_rule::IntegRule,
                              mo::MaterialOrientation)
    self=  FEMMBase{T} (fes, integration_rule, mo)
end

# Construct with the default orientation matrix updater.
function  FEMMBase{T<:FESet} (fes::T,
                    integration_rule::IntegRule)
    self=  FEMMBase{T} (fes, integration_rule, MaterialOrientation())
end

#  Associated geometry field. Default is there is none, so any
#  operation that requires a geometry field needs to be supplied it.
#  There may be operations that could benefit from pre-computations
#  that involve a geometry field. If so, associating the geometry
#  field gives the FEMM a chance to save on repeated computations.
#         assoc_geom = []; % associated geometry field
#     end

function  integrationdata (self::FEMMBase)
    # Calculate the data needed for  numerical quadrature.
    
    pc::JFFltMat = self.integration_rule.param_coords;
    w::JFFltMat  =  self.integration_rule.weights ;
    npts::JFInt = self.integration_rule.npts;
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns= Array(JFFltMat,1,npts);
    gradNparams= Array(JFFltMat,1,npts);
    for j=1:npts
        Ns[j] = FESetModule.bfun(self.fes,pc[j,:]);
        gradNparams[j] = FESetModule.bfundpar(self.fes,pc[j,:]);
    end
    return npts::JFInt, Ns::Array{JFFltMat,2}, gradNparams::Array{JFFltMat,2}, w::JFFltMat, pc::JFFltMat
end
export integrationdata

#         function fen2fe_map =node_to_element_map(self)
#      Create a map from nodes to elements.
#     
#      function fen2fe_map =node_to_element_map(self)
#     
#             if (isempty(self.fen2fe_map))
#                 self.fen2fe_map=fenode_to_fe_map (struct ('fes',self.fes));
#             end
#             fen2fe_map =self.fen2fe_map;
#         end


function integratefieldfunction{T<:Number,R} (self::FEMMBase,
                                              geom::NodalField{JFFlt},
                                              afield::NodalField{T},
                                              fh::Function,
                                              initial::R;
                                              m::JFInt=-1)
    #      Integrate a nodal-field function over the discrete manifold.
    #     
    #      function result = integrate_function (self, geom, a_field, fh, varargin)
    #     
    #      Integrate given function of a given field over the finite elements
    #      covering a manifold.
    #     
    #      Arguments:
    #      self = femm, geom = geometry field, afield = an arbitrary field,
    #      fh   = function f(Location,Field_value_at_location)
    #      m = optional: manifold dimension to be supplied to Jacobianmdim().
    
    fes=self.fes                # finite elements
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(afield); # number of degrees of freedom per node
    const nne::JFInt = nfensperfe(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self);
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    a::JFFltMat =zeros(JFFlt,nne,ndn); # array of field DOFS-- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    val::JFFltMat =zeros(JFFlt,1,ndn); # field value at the point -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    if m >= 0
        # Either the manifold dimension was supplied
    else
        m=mdim;# ...Or it is implied
    end
    result = initial;           # initial value for the result
    for i=1:nfes #Now loop over all fes in the block
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        gathervaluesasmat!(afield,a,conn);# retrieve element dofs
        for j=1:npts #Loop over all integration points
            At_mul_B!(loc,Ns[j],x);# Quadrature point location
            At_mul_B!(val,Ns[j],a);# Field value at the quadrature point
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianmdim(fes,conn, Ns[j], J, x, m);# Jacobian
            result = result + fh(loc,val)*Jac*w[j];
        end
    end
    return result
end
export integratefieldfunction

function integratefunction (self, geom, fh, varargin)
    #      Integrate a function over the discrete manifold.
    #    
    #      Integrate some scalar function over the geometric cells. When the scalar
    #      function returns just +1 [measure(femm,geom,inline('1'))], the result
    #      measures the volume (number of points, length, area, 3-D volume,
    #      according to the manifold dimension). When the function returns the mass
    #      density, the method measures the mass, when the function returns the
    #      x-coordinate equal measure the static moment with respect to the y- axis,
    #      and so on.
    #     
    #      Arguments: self = femmlock, geom = geometry field, fh   = function handle
    #      varargin= optional: manifold dimension to be supplied to Jacobian_mdim().
    #     
    #      Example:
    #          V=measure(femm,geom,@(x)(1))% Volume
    #          S=measure(femm,geom,@(x)(x))% Static moments
    #          CG=S/V% Center of gravity
    #          % Now compute the moments of inertia
    #          I=measure(femm,geom,@(x)(norm(x-CG)^2*eye(3)-(x-CG)'*(x-CG)))
    #          mass=V*rhos;
    #          Inertia=I*rhos;
    #             fes =self.fes;
    #             if nargin >=4
    #                 m =varargin{1};% Either the manifold dimension was supplied
    #             else
    #                 m= fes.dim;% ...Or it is implied
    #             end

    fes=self.fes                # finite elements
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(P); # number of degrees of freedom per node
    const nne::JFInt = nfensperfe(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self);
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    result = 0.00;# Initialize the result
    for i=1:nfes  # Now loop over all fes in the set
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        for j=1:npts #Loop over all integration points
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianmdim(fes,conn, Ns[j], J, x);# Jacobian
            result = result + fh(loc)*Jac*w(j);
        end
    end
    return result
end
export integratefunction

#         function self=associate_geometry(self,geom)
#      Associate a geometry field with the FEMM.
#     
#              function self=associate_geometry(self,geom)
#     
#      geom= Geometry field to associate with self. Pass in []
#         (empty array) to disassociate any existing geometry from the
#         self FEMM.
#     
#      Default is there is none, so any
#      operation that requires a geometry field needs to be supplied it.
#      There may be operations that could benefit from pre-computations
#      that involve a geometry field. If so, associating the geometry
#      field gives the FEMM a chance to save on repeated computations.
#             self.assoc_geom = geom;
#         end

function distribloads{S<:FESet,T<:Number}(self::FEMMBase{S}, 
                      geom::NodalField{JFFlt},
                      P::NodalField{T},
                      fi::ForceIntensity{T},
                                          m::JFInt)
    assembler = SysvecAssembler(0.0*P.values[1])#T(0.0))
    return distribloads(self,assembler,geom,P,fi,m)
end
export distribloads

function distribloads{S<:FESet,T<:Number,A<:SysvecAssemblerBase}(self::FEMMBase{S}, assembler::A,
                      geom::NodalField{JFFlt},
                      P::NodalField{T},
                      fi::ForceIntensity{T},
                      m::JFInt)
    #  % Compute the distributed-load vector.
    #  %
    #  % Return the assembled vector due to either internal heat generation
    #  % (load per unit volume), or due to applied heat flux on the surface.
    #  %
    #  % Arguments
    #  %       assembler =  descendent of sysvec_assembler
    #  %       geom=geometry field
    #  %       P=perturbation pressure field
    #  %       fi=force intensity object
    #  %       m= manifold dimension, 2= surface, 3= volume

    fes=self.fes;         # the finite element set
    # Constants
    const nfes::JFInt= count(fes); # number of finite elements in the set
    const ndn::JFInt= ndofn(P); # number of degrees of freedom per node
    const nne::JFInt = nfensperfe(fes); # number of nodes per element
    const sdim::JFInt = ndofn(geom);            # number of space dimensions
    const mdim::JFInt= manifdim(fes);     # manifold dimension of the element
    const Cedim::JFInt =ndn*nne;                          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self);
    # Prepare some data:
    conn::JFIntMat=zeros(JFInt,nne,1); # element nodes -- used as a buffer
    x::JFFltMat =zeros(JFFlt,nne,sdim); # array of node coordinates -- used as a buffer
    dofnums::JFIntMat=zeros(JFInt,1,Cedim); # degree of freedom array -- used as a buffer
    loc::JFFltMat =zeros(JFFlt,1,sdim); # quadrature point location -- used as a buffer
    J::JFFltMat =eye(JFFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
    Fe::JFMat{T} =zeros(T,Cedim,1);  
    startassembly!(assembler, P.nfreedofs);
    for i=1:nfes # Loop over elements  
        getconn!(fes,conn,i);
        gathervaluesasmat!(geom,x,conn);# retrieve element coordinates
        fill!(Fe, 0.0);
        for j=1:npts
            At_mul_B!(loc,Ns[j],x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix 
            Jac = FESetModule.Jacobianvolume(fes,conn, Ns[j], J, x);# Jacobian
            fi=updateforce!(fi,loc,J,getlabel(fes,i)); # retrieve the applied load
            Factor::JFFlt = (Jac * w[j]);
            rx::JFInt=1;
            for kx=1:nne # all the nodes
                for mx=1:ndn   # all the degrees of freedom
                    Fe[rx] = Fe[rx] + Ns[j][kx] * fi.force[mx]*Factor;
                    rx=rx+1;    # next component of the vector
                end    
            end 
        end
        gatherdofnumsasvec!(P,dofnums,conn);
        assemble!(assembler, Fe, dofnums); 
    end
    F= makevector!(assembler);
    return F
end
export distribloads

end
