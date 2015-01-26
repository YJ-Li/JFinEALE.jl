module MaterialDeformationLinearModule

using JFinEALE.JFFoundationModule
using JFinEALE.DeformationModelReductionModule
using JFinEALE.MaterialDeformationModule
using JFinEALE.PropertyDeformationLinearModule
using JFinEALE.PropertyDeformationLinearIsoModule
using JFinEALE.PropertyDeformationLinearOrthoModule

# % Class of linear elastic material for a variety of model reduction strain/stress states.
# %
# % This class represents deformable materials, linear elastic
# strain/stress states.           # 

type  MaterialDeformationLinear{P<:PropertyDeformationLinear}
    property::P
end
export MaterialDeformationLinear


################################################################################
# 3-D solid model

function tangentmoduli!{P<:PropertyDeformationLinear,
    MR<:DeformationModelReduction3D}(::Type{MR},
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
    tangentmoduli3d!(self.property, D; context...);
    return D
end
export tangentmoduli!

function update!{MR<:DeformationModelReduction3D}(::Type{MR},
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
    D=zeros(JFFlt,6,6)
    tangentmoduli!(MR, self, D; context...);# local material stiffness
    tSigma = thermalstress(MR, self; context...);# stress in local coordinates
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

function thermalstress{MR<:DeformationModelReduction3D}(::Type{MR},
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
            if (dT!= 0.0)    # nonzero temperature differential?
                D=zeros(JFFlt,6,6)
                tangentmoduli!(MR, self, D; context...);# local material stiffness
                v = -D[:,1:3]*(dT*self.property.CTE);
                return v
            end
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
    
    D3d=zeros(JFFlt,6,6)
    tangentmoduli3d!(self.property, D3d; context...);
    const ix=[1, 2, 4];
    for i=1:length(ix)
        for j=1:length(ix)
            D[j,i]=D3d[ix[j],ix[i]];
        end        
    end
    return D
end
export tangentmoduli!


function update!{MR<:DeformationModelReduction2DStrain}(::Type{MR},
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
        #           outputRm=optional orientation matrix in which output should 
        #               supplied   
        #
        #   It is assumed that stress is output in m-component vector 
        #           form. m=3 for plane stress, m=4 for plane strain or axially 
        #           symmetric
        #   The output arguments are
        #     out=requested quantity
        #     newms=new material state; don't forget that if the update is final
        #           the material state newms must be assigned and stored.  Otherwise
        #           the material update is lost!
    #
    Ev=JFFlt[]                  # it is empty: we need to get it from context
    output=:Cauchy
    dT= 0.0
    for arg in context
        sy, val = arg
        if sy==:strain
            Ev=val
        elseif sy==:output
            output=val
        elseif sy==:dT
            dT=val
        end
    end
    D=zeros(JFFlt,3,3)
    tangentmoduli!(MR,self,D; context...)
    tSigma = thermalstress(MR,self; context...);
    stress = D * Ev + tSigma;
    # Need to output the along-the-thickness  stress
    CTE=self.property.CTE
    D3d=zeros(JFFlt,6,6)
    tangentmoduli3d!(self.property, D3d; context...);
    sz=D3d[3,1:2]*Ev[1:2]-dT[1]*D3d[3,1:2]*CTE[1:2]-dT[1]*D3d[3,3]*CTE[3];
    # sigmax, sigmay, tauxy, sigmaz
    stress = [vec(stress[1:3]), vec(sz)];
    if output==:Cauchy
        out = stress;
    elseif output==:pressure
        out = -sum(stress[[1,2,4]])/3;
    elseif output==:volstrain
        out = sum(Ev[1:2]);  
    elseif output==:princCauchy
        t=zeros(JFFlt,3,3)
        t = stress4vto3x3t!(stress,t);
        ep=eig(t);
        out =sort(ep[1]);    
    else
        out = stress;
    end
    newms = ms;
    return out,newms;
end

function thermalstress{MR<:DeformationModelReduction2DStrain}(::Type{MR},
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
            D=zeros(JFFlt,3,3)
            tangentmoduli!(MR, self, D; context...);# local material stiffness
            # there is no sigmaz in this vector!
            CTE=self.property.CTE
            v = -dT[1]*[D[1:2,1:2]*CTE[1:2]+CTE[3]*D[1:2,3]; 0];
            return v
        end
    end
    return  zeros(JFFlt, 3, 1);
end
export thermalstress

################################################################################
# Plane stress model

function tangentmoduli!{P<:PropertyDeformationLinear,
    MR<:DeformationModelReduction2DStress}(::Type{MR},
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

    D3d=zeros(JFFlt,6,6)
    tangentmoduli3d!(self.property, D3d; context...);
    Dt = D3d[1:2, 1:2]- D3d[1:2,3]*D3d[3,1:2]/D3d[3,3];
    const ix=[1, 2, 4];
    for i=1:2
        for j=1:2
            D[j,i]= Dt[j,i];
        end        
    end
    for i=1:3
        D[3,i]= D[i,3]= D3d[4,ix[i]];
    end
    return D
end
export tangentmoduli!

function update!{MR<:DeformationModelReduction2DStress}(::Type{MR},
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
        #           outputRm=optional orientation matrix in which output should 
        #               supplied   
        #
        #   It is assumed that stress is output in m-component vector 
        #           form. m=3 for plane stress, m=4 for plane strain or axially 
        #           symmetric
        #   The output arguments are
        #     out=requested quantity
        #     newms=new material state; don't forget that if the update is final
        #           the material state newms must be assigned and stored.  Otherwise
        #           the material update is lost!
    #
    Ev=JFFlt[]                  # it is empty: we need to get it from context
    output=:Cauchy
    for arg in context
        sy, val = arg
        if sy==:strain
            Ev=val
        elseif sy==:output
            output=val
        end
    end
    D=zeros(JFFlt,3,3)
    tangentmoduli!(MR,self,D; context...)
    tSigma = thermalstress(MR,self; context...);
    stress = D * Ev + tSigma;
    if output==:Cauchy
        out = stress;
    elseif output==:pressure
        out = -sum(stress(1:2))/3;
    elseif output==:volstrain
        out = sum(Ev(1:2));     # actually this is probably incorrect
        #for plane stress: the transverse strain is not accounted for
    elseif output==:princCauchy
        t=zeros(JFFlt,3,3)
        t = stress3vto3x3t!(stress,t);
        ep=eig(t);
        out =sort(ep[1]);    
    else
        out = stress;
    end
    newms = ms;
    return out,newms;
end

function thermalstress{MR<:DeformationModelReduction2DStress}(::Type{MR},
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
            D=zeros(JFFlt,3,3)
            tangentmoduli!(MR, self, D; context...);# local material stiffness
            v = -D*(dT[1]*[self.property.CTE[1:2],0.0]);
            return v
            # switch self.reduction
            #     case 'axisymm'
            #         D = tangent_moduli(self, context);
            #         v = -D*context.dT*[alphas(1:3).*ones(3, 1); 0];
            #     case 'strain'
            #         D=  self.property.tangent_moduli(context);% need 3-D
            #         v = -context.dT*[D(1:2, 1:2)*(alphas(1:2).*ones(2,1))+...
            #             alphas(3)*D(1:2,3); 0];
            
        end
    end
    return  zeros(JFFlt, 3, 1);
end
export thermalstress

            

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
    
    D3d=zeros(JFFlt,6,6)
    tangentmoduli3d!(self.property, D3d; context...);
    for i=1:4
        for j=1:4
            D[j,i]= D3d[j,i];
        end        
    end
    return D
end
export tangentmoduli!

function update!{MR<:DeformationModelReduction2DAxisymm}(::Type{MR},
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
        #           outputRm=optional orientation matrix in which output should 
        #               supplied   
        #
        #   It is assumed that stress is output in m-component vector 
        #           form. m=3 for plane stress, m=4 for plane strain or axially 
        #           symmetric
        #   The output arguments are
        #     out=requested quantity
        #     newms=new material state; don't forget that if the update is final
        #           the material state newms must be assigned and stored.  Otherwise
        #           the material update is lost!
    #
    Ev=JFFlt[]                  # it is empty: we need to get it from context
    output=:Cauchy
    for arg in context
        sy, val = arg
        if sy==:strain
            Ev=val
        elseif sy==:output
            output=val
        end
    end
    D=zeros(JFFlt,4,4)
    tangentmoduli!(MR,self,D; context...)
    tSigma = thermalstress(MR,self; context...);
    stress = D * Ev + tSigma;
    # println("D=$( D  )")
    #  println("Ev=$( Ev )")
    # println("stress=$( stress )")
    if output==:Cauchy
        out = stress;
    elseif output==:pressure
        out = -sum(stress[[1,2,4]])/3;
    elseif output==:volstrain
        out = sum(Ev[[1,2,4]]);  
    elseif output==:princCauchy
        t=zeros(JFFlt,3,3)
        t = stress4vto3x3t!(stress[[1,2,4,3]],t);
        ep=eig(t);
        out =sort(ep[1]);    
    else
        out = stress;
    end
    newms = ms;
    return out,newms;
end

function thermalstress{MR<:DeformationModelReduction2DAxisymm}(::Type{MR},
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
            D=zeros(JFFlt,4,4)
            tangentmoduli!(MR, self, D; context...);# local material stiffness
            v = -D*(dT[1]*[self.property.CTE; 0.0]);
            return v
            #     case 'strain'
            #         D=  self.property.tangent_moduli(context);% need 3-D
            #         v = -context.dT*[D(1:2, 1:2)*(alphas(1:2).*ones(2,1))+...
            #             alphas(3)*D(1:2,3); 0];
            
        end
    end
    return  zeros(JFFlt, 4, 1);
end
export thermalstress

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
    
    D3d=zeros(JFFlt,6,6)
    tangentmoduli3d!(self.property, D3d; context...);
    D[1,1] = D3d[1, 1]- D3d[1,2:3]*D3d[2:3,2:3]\D3d[2:3,1];
    return D
end
export tangentmoduli!

end
