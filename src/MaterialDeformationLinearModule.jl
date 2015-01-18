module MaterialDeformationLinearModule

using JFinEALE.JFFoundationModule
using JFinEALE.DeformationModelReductionModule
using JFinEALE.MaterialDeformationModule
using JFinEALE.PropertyDeformationLinearModule
using JFinEALE.PropertyDeformationLinearIsoModule

# % Class of linear elastic material for a variety of model reduction strain/stress states.
# %
# % This class represents deformable materials, linear elastic
# strain/stress states.           # 

type  MaterialDeformationLinear{P<:PropertyDeformationLinear}
    property::P
    bufferedD::JFFltMat
end
export MaterialDeformationLinear

function MaterialDeformationLinear{P<:PropertyDeformationLinear}(property::P)
    MaterialDeformationLinear(property,zeros(JFFlt,6,6));
end


################################################################################
# 3-D solid model

function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeformationModelReduction3D},
                        self::MaterialDeformationLinear{P}, 
                        D::JFFltMat; context...)
    # # Calculate the material stiffness matrix.
    # # 
    # # Arguments
    # #     m=material
    # #     context=structure with mandatory and optional fields that are required 
    # # by the specific material implementation.
    # #
    # # the output arguments are
    # #     D=matrix 6x6 in the local material orientation Cartesian basis
    # #
    # #
    tangentmoduli!(self.property, D; context...);
end
export tangentmoduli!

function update! (::Type{DeformationModelReduction3D},
                        self::MaterialDeformationLinear, ms; context...)
# Update material state.
#
# function [out, newms] = update (self, ms, context)
#
# Update material state.  Return the updated material state, and the
# requested quantity (default is the stress).
#   Call as:
#     [out,newms] = update(m, ms, context)
#  where
#     m=material
#     ms = material state
#     context=structure
#        with mandatory fields
#           strain=strain vector  in the local material
#               directions (which may be the same as the global coordinate
#               directions)
#        and optional fields
#           output=type of quantity to output, and interpreted by the
#               particular material; [] is returned when the material does not
#               recognize the requested quantity to indicate uninitialized
#               value.  It can be tested with isempty ().
#                  output ='Cauchy' - Cauchy stress; this is the default
#                      when output type is not specified.
#                  output ='princCauchy' - principal Cauchy stress;
#                  output ='pressure' - pressure;
#                  output ='vol_strain' - volumetric strain;
#                  output ='vonMises' - von Mises stress;  
#           outputRm=optional orientation matrix in which output should 
#               supplied   
#
#   It is assumed that stress is output in 6-component vector form.
#   The output arguments are
#     out=requested quantity
#           Remember: the output is expressed in the local material 
#           orientation  coordinates.
#     newms=new material state; don't forget that if the update is final
#           the material state newms must be assigned and stored.  Otherwise
#           the material update is lost!
    #
    
# if (isfield(context,'strain'))
#     Ev = context.strain;# strain in material coordinates
# else# This is an approximation valid only for small displacements
#     gradu=context.F-eye(3);
#     Ev = strain_3x3t_to_6v (self,(gradu+gradu')/2);
# end
    tangentmoduli!(self, self.bufferedD; context...);# local material stiffness
    tSigma = thermalstress(self; context...);# stress in local coordinates
    stress = D * Ev + tSigma;
    # if isfield(context,'output')
    #     switch context.output
    #         case 'Cauchy'
    #             out = stress;
    #         case 'vol_strain'
    #             out = sum(Ev(1:3));
    #         case 'pressure'
    #             out = -(sum(stress(1:3))/3);
    #         case 'princCauchy'
    #             t = stress_6v_to_3x3t (self,stress);
    #             [V,D]=eig(t);
    #             out =sort(diag(D),'descend');
    #         case {'vonMises','vonmises','von_mises','vm'}
    #             s1=stress(1);s2=stress(2);s3=stress(3);
    #             s4=stress(4);s5=stress(5);s6=stress(6);
    #             out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
    #         otherwise
    #             out = [];
    #     end
    # else
    #     out = stress;
    # end
    newms = ms;
    return;
end
export update!

function thermalstress(::Type{DeformationModelReduction3D},
                        self::MaterialDeformationLinear; context...)
# Calculate vector of thermal stress components.
#
# function v = thermal_stress(self,context)
#
#   Call as:
#     v = thermal_stress(m,context)
#  where
#     m=material
#     context=structure; see the update() method
    #
    for arg in args
        sy, val = arg
        if sy==:dT
            dT=val
            tangentmoduli!(self, self.bufferedD; context...);# local material stiffness
            v = -self.bufferedD[:,1:3]*(dT*CTE);
            return v
        end
    end
    return  zeros(JFFlt, 6, 1);
end
export thermalstress

################################################################################
# Plane strain model

function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeformationModelReduction2DStrain},
                        self::MaterialDeformationLinear{P}, 
                        D::JFFltMat; context...)
    # # Calculate the material stiffness matrix.
    # # 
    # # Arguments
    # #     m=material
    # #     context=structure with mandatory and optional fields that are required 
    # # by the specific material implementation.
    # #
    # # the output arguments are
    # #     D=matrix 6x6 in the local material orientation Cartesian basis
    # #
    
    
    tangentmoduli!(self.property, self.bufferedD; context...);
    const ix=[1, 2, 4];
    for i=1:3
        for j=1:3
            D[j,i]=self.bufferedD[ix[j],ix[i]];
        end        
    end    
end
export tangentmoduli!


function thermalstress(::Type{DeformationModelReduction2DStrain},
                        self::MaterialDeformationLinear; context...)
# Calculate vector of thermal stress components.
#
# function v = thermal_stress(self,context)
#
#   Call as:
#     v = thermal_stress(m,context)
#  where
#     m=material
#     context=structure; see the update() method
    #
    for arg in context
        sy, val = arg
        if sy==:dT
            dT=val
            tangentmoduli!(self, self.bufferedD; context...);# local material stiffness
            v = -self.bufferedD[:,1:3]*(dT*CTE);
            return v
        end
    end
    return  zeros(JFFlt, 6, 1);
end
export thermalstress

################################################################################
# Plane stress model

function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeformationModelReduction2DStress},
                        self::MaterialDeformationLinear{P}, 
                        D::JFFltMat; context...)
    # # Calculate the material stiffness matrix.
    # # 
    # # Arguments
    # #     m=material
    # #     context=structure with mandatory and optional fields that are required 
    # # by the specific material implementation.
    # #
    # # the output arguments are
    # #     D=matrix 6x6 in the local material orientation Cartesian basis
    # #
    
    tangentmoduli3d!(self.property, self.bufferedD; context...);
    Dt = self.bufferedD[1:2, 1:2]- self.bufferedD[1:2,3]*self.bufferedD[3,1:2]/self.bufferedD[3,3];
    const ix=[1, 2, 4];
    for i=1:2
        for j=1:2
            D[j,i]= Dt[j,i];
        end        
    end
    for i=1:3
        D[3,i]= D[i,3]= self.bufferedD[4,ix[i]];
    end    
end
export tangentmoduli!

################################################################################
# Axially symmetric model

function tangentmoduli!(::Type{DeformationModelReduction2DAxisymm},
                        self::MaterialDeformationLinear, 
                        D::JFFltMat; context...)
    # # Calculate the material stiffness matrix.
    # # 
    # # Arguments
    # #     m=material
    # #     context=structure with mandatory and optional fields that are required 
    # # by the specific material implementation.
    # #
    # # the output arguments are
    # #     D=matrix 6x6 in the local material orientation Cartesian basis
    # #
    
    
    tangentmoduli3d!(self.property, self.bufferedD; context...);
    for i=1:4
        for j=1:4
            D[j,i]= self.bufferedD[j,i];
        end        
    end
end
export tangentmoduli!

################################################################################
# 1D model

function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeformationModelReduction1D},
                        self::MaterialDeformationLinear{P}, 
                        D::JFFltMat; context...)
    # # Calculate the material stiffness matrix.
    # # 
    # # Arguments
    # #     m=material
    # #     context=structure with mandatory and optional fields that are required 
    # # by the specific material implementation.
    # #
    # # the output arguments are
    # #     D=matrix 6x6 in the local material orientation Cartesian basis
    # #
    
    
    tangentmoduli3d!(self.property, self.bufferedD; context...);
    D[1,1] = self.bufferedD[1, 1]- self.bufferedD[1,2:3]*self.bufferedD[2:3,2:3]\self.bufferedD[2:3,1];
end
export tangentmoduli!

end
