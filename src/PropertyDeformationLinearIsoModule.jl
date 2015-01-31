module PropertyDeformationLinearIsoModule


using JFinEALE.JFFoundationModule
using JFinEALE.PropertyDeformationLinearModule
import JFinEALE.PropertyDeformationLinearModule.tangentmoduli3d!


# Class of properties for linearly elastic isotropic homogeneous materials. 

type PropertyDeformationLinearIso<:PropertyDeformationLinear
    mass_density::JFFlt                  # mass density
    E::JFFlt
    nu::JFFlt
    CTE::JFFltVec               # three thermal expansion coefficients
    D::JFFltMat
end
export PropertyDeformationLinearIso

const mI = diagm([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
const m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];

function PropertyDeformationLinearIso(E::JFFlt,nu::JFFlt)
    rho= 0.0;
    if (E<=0.0)
        error("Non-positive Young's modulus!");
    end
    if (nu<0.0)
        error("Negative Poisson ratio!");
    end
    lambda = E * nu / (1 + nu) / (1 - 2*(nu));
    mu = E / 2. / (1+nu);
    D = lambda * m1 * m1' + 2. * mu * mI;
    return PropertyDeformationLinearIso(rho,E,nu,zeros(JFFlt,3),D)
end


function PropertyDeformationLinearIso(rho::JFFlt,E::JFFlt,nu::JFFlt)
    if (rho<=0.0)
        error("Non-positive mass density!");
    end
    if (E<=0.0)
        error("Non-positive Young's modulus!");
    end
    if (nu<0.0)
        error("Negative Poisson ratio!");
    end
    lambda = E * nu / (1 + nu) / (1 - 2*(nu));
    mu = E / 2. / (1+nu);
    D = lambda * m1 * m1' + 2. * mu * mI;
    return PropertyDeformationLinearIso(rho,E,nu,zeros(JFFlt,3),D)
end

function PropertyDeformationLinearIso(rho::JFFlt,E::JFFlt,nu::JFFlt,CTE::JFFlt)
    if (rho<=0.0)
        error("Non-positive mass density!");
    end
    if (E<=0.0)
        error("Non-positive Young's modulus!");
    end
    if (nu<0.0)
        error("Negative Poisson ratio!");
    end
    lambda = E * nu / (1 + nu) / (1 - 2*(nu));
    mu = E / 2. / (1+nu);
    D = lambda * m1 * m1' + 2. * mu * mI;
    return PropertyDeformationLinearIso(rho,E,nu,zeros(JFFlt,3)+CTE,D)
end

function PropertyDeformationLinearIso(;rho::JFFlt= 1.0,E::JFFlt= 1.0,nu::JFFlt= 0.0,CTE::JFFlt=  0.0)
    return PropertyDeformationLinearIso(rho,E,nu,CTE)
end

function  tangentmoduli3d!(self::PropertyDeformationLinearIso,D::JFFltMat; context...)
    # Calculate the material stiffness matrix.
    #
    # function val = tangent_moduli(self, context)
    #
    # Arguments
    #    self=property
    #    context=structure with optional field:
    #                kind = string with possible values 'lambda',
    #                       'lambda_shear', 'bulk', 'bulk_shear'
    #
    # Output
    #     D=matrix 6x6 in the local material orientation Cartesian basis

    # kind=:full;
    # for arg in context
    #     sy, val = arg
    #     if sy==:kind
    #         dT=val
    #         tangentmoduli!(self, self.bufferedD; context...);# local material stiffness
    #         v = -self.bufferedD[:,1:3]*(dT*CTE);
    #         return v
    #     end
    # end
    # if kind!=:full
    #      error("Unknown kind")
    #     # switch context.kind
    #     # case 'lambda'
    #     # val = tangent_moduli_lambda(self, context);
    #     # case 'lambda_shear'
    #     # val = tangent_moduli_lambda_shear(self, context);
    #     # case 'bulk'
    #     # val = tangent_moduli_bulk(self, context);
    #     # case 'bulk_shear'
    #     # val = tangent_moduli_bulk_shear(self, context);
    #     # case 'constrained'
    #     # #val = tangent_moduli_bulk(self, context);
    #     # val = tangent_moduli_lambda(self, context);
    #     # case 'unconstrained'
    #     # #val = tangent_moduli_bulk_shear(self, context);
    #     # val = tangent_moduli_lambda_shear(self, context);
    #     # otherwise
    #     # val = [];
    # end
    #println("in tangentmoduli3d!: $(self.D)")
    copy!(D,self.D);
end
export   tangentmoduli3d!

end
    
    
#     methods (Hidden, Access = private)
        
#         function D = tangent_moduli_lambda(self, context)
#         # Calculate the part of the material stiffness matrix that corresponds to
#         # the lambda Lame coefficient.
#             E = self.E;
#             nu = self.nu;
#             lambda = E * nu / (1 + nu) / (1 - 2*(nu));
#             m1 = [1 1 1 0 0 0]';
#             D = lambda * m1 * m1';
#             return;
#         end
        
#         function D = tangent_moduli_lambda_shear(self, context)
#         # Calculate the part of the material stiffness matrix that correspond to shear.
#         # Note: makes sense only for isotropic materials.
#             E = self.E;
#             nu = self.nu;
#             mu     = E / (2 * (1 + nu));
#             mI = diag([1 1 1 0.5 0.5 0.5]);
#             D = 2 * mu * mI;
#             return;
#         end
        
#         function D = tangent_moduli_bulk(self, context)
#         # Calculate the part of the material stiffness matrix that corresponds to
#         # the bulk modulus.
#             E = self.E;
#             nu = self.nu;
#             B = E / 3 / (1 - 2*(nu));
#             m1 = [1 1 1 0 0 0]';
#             D = B * m1 * m1';
#             return;
#         end
        
#         function D = tangent_moduli_bulk_shear(self, context)
#         # Calculate the part of the material stiffness matrix that correspond to shear.
#         # Note: makes sense only for isotropic materials.
#             E = self.E;
#             nu = self.nu;
#             mu     = E / (2 * (1 + nu));
#             D = mu * [2/3*[2 -1 -1; -1 2 -1; -1 -1 2] zeros(3,3); zeros(3,3) eye(3,3) ];
#             return;
#         end
        
#     end
    
#     # Some checking:
#     #     mI = diag([1 1 1 0.5 0.5 0.5]);
#     # m1 = [1 1 1 0 0 0]';
#     # mId =eye(6) - 1/3*m1 * m1';
#     # K  = sym('K','real');
#     # G  = sym('G','real');
#     # D = K*m1 * m1'+2*G*mI*mId
#     #  (1/3*m1 * m1')*D*(1/3*m1 * m1')
#     #  mId*D*mId
#     #  (1/3*m1 * m1')*D*(1/3*m1 * m1')+ mId*D*mId-D
# end


# end
