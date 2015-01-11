module PropertyAcousticFluidModule


using JFFoundationModule

# Type for acoustic fluid material property.
# Fields:
#     bulk_modulus=bulk modulus
#     mass_density=mass density
#
type PropertyAcousticFluid
    bulk_modulus::JFFlt;# Bulk modulus
    mass_density::JFFlt;# Mass density
end
export PropertyAcousticFluid

end
