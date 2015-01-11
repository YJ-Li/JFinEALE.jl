module ForceIntensityModule

using JFFoundationModule

type ForceIntensity{T<:Number}
    #   % Distributed force (force intensity) class.
    #     %
    #     % The force intensity class. The physical units are
    #     % force per unit volume, where volume depends on to which manifold
    #     % the force is applied:
    #     % force/length^3 (when applied to a 3-D solid),
    #     % force/length^2 (when applied to a surface),
    #     % force/length^1 (when applied along a curve),
    #     % or force/length^0 (when applied at a point).
    #    % Magnitude of the distributed force.
    #         % Either a constant vector or a function handle.
    force::JFMat{T}    # buffer where the current value force can be retrieved
    callback::Function # client's function to return the force intensity magnitude
    callbackactive::Bool
end
export ForceIntensity

function  ForceIntensity{T<:Number}(force::JFMat{T})
    #             magn=vector of magnitudes, dimension must correspond to the
    #             %            dimension of the displacement/geometry fields; magn may be an array
    #             %            of doubles (representing a constant, or uniform, force intensity)
    #             % or a function handle/in-line function with signature
    #             %                     function val=f(xyz)
    #             %            or
    #             %                     function val=f(xyz,J)
    #             %            where
    #             %                     xyz = location, at which the force intensity
    #             %                           is to be evaluated.
    #             %                     J = Jacobian matrix at the location above.
    #             %                        The Jacobian matrix could be useful for instance
    #             %                         for the calculation of the normal to the surface.

    function nullcallback(xyz,J)
        error("Callback not active");
    end
    
    return ForceIntensity(copy(force), nullcallback, false);
end

function  ForceIntensity{T<:Number}(force::JFVec{T})
    #             magn=vector of magnitudes, dimension must correspond to the
    #             %            dimension of the displacement/geometry fields; magn may be an array
    #             %            of doubles (representing a constant, or uniform, force intensity)
    #             % or a function handle/in-line function with signature
    #             %                     function val=f(xyz)
    #             %            or
    #             %                     function val=f(xyz,J)
    #             %            where
    #             %                     xyz = location, at which the force intensity
    #             %                           is to be evaluated.
    #             %                     J = Jacobian matrix at the location above.
    #             %                        The Jacobian matrix could be useful for instance
    #             %                         for the calculation of the normal to the surface.

    function nullcallback(xyz,J)
        error("Callback not active");
    end
    
    return ForceIntensity(copy(reshape(force,length(force),1)), nullcallback, false);
end

function  ForceIntensity{T<:Number}(force::T)
    #             magn=vector of magnitudes, dimension must correspond to the
    #             %            dimension of the displacement/geometry fields; magn may be an array
    #             %            of doubles (representing a constant, or uniform, force intensity)
    #             % or a function handle/in-line function with signature
    #             %                     function val=f(xyz)
    #             %            or
    #             %                     function val=f(xyz,J)
    #             %            where
    #             %                     xyz = location, at which the force intensity
    #             %                           is to be evaluated.
    #             %                     J = Jacobian matrix at the location above.
    #             %                        The Jacobian matrix could be useful for instance
    #             %                         for the calculation of the normal to the surface.

    function nullcallback(xyz,J)
        error("Callback not active");
    end
    
    return ForceIntensity(copy(reshape([force],1,1)), nullcallback, false);
end

function  ForceIntensity{T<:Number}(force::JFMat{T},callback::Function)
    #             magn=vector of magnitudes, dimension must correspond to the
    #             %            dimension of the displacement/geometry fields; magn may be an array
    #             %            of doubles (representing a constant, or uniform, force intensity)
    #             % or a function handle/in-line function with signature
    #             %                     function val=f(xyz)
    #             %            or
    #             %                     function val=f(xyz,J)
    #             %            where
    #             %                     xyz = location, at which the force intensity
    #             %                           is to be evaluated.
    #             %                     J = Jacobian matrix at the location above.
    #             %                        The Jacobian matrix could be useful for instance
    #             %                         for the calculation of the normal to the surface.

    return ForceIntensity(copy(force), callback, true);
end

function getforce!{T<:Number}(self::ForceIntensity{T},xyz::JFFltMat,J::JFFltMat)
    if (self.callbackactive)
        copy!(self.force,self.callback(xyz,J));
    else
        # be nothing to be done: the data is already in the buffer
    end
    return self
end
export getforce!

end
