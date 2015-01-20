module ForceIntensityModule

using JFinEALE.JFFoundationModule

#    Distributed force (force intensity) class.
    #     
    #      The force intensity class. The physical units are
    #      force per unit volume, where volume depends on to which manifold
    #      the force is applied:
    #      force/length^3 (when applied to a 3-D solid),
    #      force/length^2 (when applied to a surface),
    #      force/length^1 (when applied along a curve),
    #      or force/length^0 (when applied at a point).
    #     Magnitude of the distributed force.
    #          Either a constant vector or a function handle.

type ForceIntensity{T<:Number}
    g!::Function # function to update the force
    force::JFMat{T}    # buffer where the current value of the force can be retrieved
end
export ForceIntensity

# Construct force intensity when the function to compute the intensity matrix is given.
function ForceIntensity{T<:Number}(::Type{T},computeforce::Function)
    force= Array(T,0,0);       # Empty buffer, needs to be allocated on the first call
    function g!(forceout::JFMat{T}, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        r=computeforce(XYZ, tangents, fe_label)::JFMat{T}
        if (length(forceout)== 0) # If the buffer is not initialized yet
            forceout = zeros(T,size(r)); 
        end
        copy!(forceout,r)
        return forceout
    end
    return ForceIntensity(g!,force);
end

# Construct force intensity when the intensity is given as a matrix.
function ForceIntensity{T<:Number}(force::JFMat{T})
    iforce= deepcopy(force);          # fill the buffer with the given matrix
    function g!(forceout::JFMat{T}, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        return forceout            # nothing to be done here, the matrix is already in the buffer
    end
    return ForceIntensity(g!,iforce);
end

# Construct force intensity when the intensity is given as a vector.
function ForceIntensity{T<:Number}(force::JFVec{T})
    iforce= reshape(deepcopy(force),length(force),1);   # fill the buffer with the given matrix
    function g!(forceout::JFMat{T}, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        return forceout            # nothing to be done here, the matrix is already in the buffer
    end
    return ForceIntensity(g!,iforce);
end

# Construct force intensity when the intensity is given as a scalar.
function  ForceIntensity{T<:Number}(force::T)
    iforce= reshape(deepcopy([force]),1,1);   # fill the buffer with the given matrix
    function g!(forceout::JFMat{T}, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        return forceout            # nothing to be done here, the matrix is already in the buffer
    end
    return ForceIntensity(g!,iforce);
end 

# Update the force intensity matrix.
# After the return it can be retrieved from the buffer (self.force).
function updateforce!{T<:Number}(self::ForceIntensity{T},XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
    self.force=self.g!(self.force, XYZ,tangents,fe_label)
    return self
end
export updateforce!

end
