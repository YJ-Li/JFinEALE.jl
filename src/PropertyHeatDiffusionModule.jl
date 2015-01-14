module PropertyHeatDiffusionModule


using JFinEALE.JFFoundationModule

# Type for heat diffusion material property.
# Fields:
#     thermal_conductivity=matrix of conductivities; 2x2 for 2-D domains, 3x3
#                  for 3-D domains; the components are given in the local
#                  coordinate frame (femm decides what that frame is by reference to the material orientation matrix).
#     specific_heat=specific heat per unit volume
#
type PropertyHeatDiffusion
    thermal_conductivity::JFFltMat;# Thermal conductivity
    specific_heat::JFFlt;# Specific heat per unit volume
end
export PropertyHeatDiffusion

function  PropertyHeatDiffusion (thermal_conductivity,specific_heat=0.0)
    return  PropertyHeatDiffusion (thermal_conductivity,specific_heat)
end

end
