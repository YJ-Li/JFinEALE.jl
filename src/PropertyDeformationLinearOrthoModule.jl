module PropertyDeformationLinearOrthoModule


using JFinEALE.JFFoundationModule
using JFinEALE.PropertyDeformationLinearModule
import JFinEALE.PropertyDeformationLinearModule.tangentmoduli3d!


# Class of properties for linearly elastic isotropic homogeneous materials. 

type PropertyDeformationLinearOrtho<:PropertyDeformationLinear
    rho::JFFlt                  # mass density
    E1::JFFlt
    E2::JFFlt
    E3::JFFlt
    nu12::JFFlt
    nu13::JFFlt
    nu23::JFFlt
    G12::JFFlt
    G13::JFFlt
    G23::JFFlt
    CTE::JFFltVec               # three thermal expansion coefficients
    D::JFFltMat
end
export PropertyDeformationLinearOrtho

const mI = diagm([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
const m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];

function PropertyDeformationLinearOrtho(rho::JFFlt,
                                        E1::JFFlt,E2::JFFlt,E3::JFFlt,
                                        nu12::JFFlt,nu13::JFFlt,nu23::JFFlt,
                                        G12::JFFlt,G13::JFFlt,G23::JFFlt,
                                        CTE1::JFFlt,CTE2::JFFlt,CTE3::JFFlt)
    if (E1<=0.0)
        error("Non-positive Young's modulus E1!");
    end
    if (E2<=0.0)
        error("Non-positive Young's modulus E2!");
    end
    if (E3<=0.0)
        error("Non-positive Young's modulus E3!");
    end
    if (G12<=0.0)
        error("Non-positive shear modulus G12!");
    end
    if (G13<=0.0)
        error("Non-positive shear modulus G13!");
    end
    if (G23<=0.0)
        error("Non-positive shear modulus G23!");
    end
    C = [  1.0/E1      -nu12/E1    -nu13/E1  0.0   0.0   0.0;
         -nu12/E1      1.0/E2      -nu23/E2  0.0   0.0   0.0;
         -nu13/E1    -nu23/E2       1.0/E3   0.0   0.0   0.0;
         0.0           0.0           0.0   1/G12   0.0   0.0;
         0.0           0.0           0.0     0.0 1/G13   0.0;
         0.0           0.0           0.0     0.0   0.0  1/G23];
    D = inv(C)
    if (rank(D)<6)
        error("Non-positive definite D!");
    end
    return PropertyDeformationLinearOrtho(rho,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,[CTE1,CTE2,CTE3],D)
end

function PropertyDeformationLinearOrtho(rho::JFFlt,
                                        E1::JFFlt,E2::JFFlt,E3::JFFlt,
                                        nu12::JFFlt,nu13::JFFlt,nu23::JFFlt,
                                        G12::JFFlt,G13::JFFlt,G23::JFFlt)
    CTE1=CTE2=CTE3= 0.0         # default
    return PropertyDeformationLinearOrtho(rho,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,CTE1,CTE2,CTE3)
end

function PropertyDeformationLinearOrtho(E1::JFFlt,E2::JFFlt,E3::JFFlt,
                                        nu12::JFFlt,nu13::JFFlt,nu23::JFFlt,
                                        G12::JFFlt,G13::JFFlt,G23::JFFlt)
    rho= 0.0                    # default
    return PropertyDeformationLinearOrtho(rho,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
end

function PropertyDeformationLinearOrtho(E1::JFFlt,E2::JFFlt,E3::JFFlt,
                                        nu12::JFFlt,nu13::JFFlt,nu23::JFFlt,
                                        G12::JFFlt,G13::JFFlt,G23::JFFlt,
                                        CTE1::JFFlt,CTE2::JFFlt,CTE3::JFFlt)
    rho= 0.0                    # default
    return PropertyDeformationLinearOrtho(rho,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,CTE1,CTE2,CTE3)
end

function  tangentmoduli3d!(self::PropertyDeformationLinearOrtho,D::JFFltMat; context...)
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

    copy!(D,self.D);
end
export   tangentmoduli3d!

end
    
   
