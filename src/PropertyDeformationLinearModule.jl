module PropertyDeformationLinearModule


using JFinEALE.JFFoundationModule


# Class of properties for linearly elastic isotropic materials. 

abstract PropertyDeformationLinear
export PropertyDeformationLinear


function  tangentmoduli3d!(self::PropertyDeformationLinear,D::JFFltMat; context...)
    error("Shouldn't have gotten here!")
end
export tangentmoduli3d!

end
