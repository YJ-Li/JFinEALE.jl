module MaterialOrientationModule

using JFinEALE.JFFoundationModule

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

type MaterialOrientation
    g!::Function # function to update the material orientation matrix
    Rm::JFFltMat # the material orientation matrix (buffer)
end
export MaterialOrientation

# Construct material orientation when the function to compute the rotation matrix is given.
function MaterialOrientation(computeRm::Function)
    Rm= Array(JFFlt,0,0);       # Empty buffer, needs to be allocated on the first call
    function g!(Rmout::JFFltMat, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        r=computeRm(XYZ, tangents, fe_label)::JFFltMat
        if (length(Rmout)== 0) # If the buffer is not initialized yet
            Rmout = zeros(size(r)); 
        end
        copy!(Rmout,r)
        return Rmout
    end
    return MaterialOrientation(g!,Rm);
end

# Construct material orientation when the rotation matrix is given.
function MaterialOrientation(Rm::JFFltMat)
    iRm= deepcopy(Rm);          # fill the buffer with the given matrix
    function g!(Rmout::JFFltMat, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        return Rmout            # nothing to be done here, the matrix is already in the buffer
    end
    return MaterialOrientation(g!,iRm);
end

# Construct default material orientation (identity matrix).
function MaterialOrientation()
    Rm= Array(JFFlt,0,0);       # Empty buffer, needs to be allocated on the first call
    function g!(Rmout::JFFltMat, XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
        if (length(Rmout)== 0) # If the buffer is not initialized yet
            Rmout = zeros(size(XYZ,2),size(tangents,2)); # space dimension x manifold dimension
            fill!(Rmout,0.0)                             # build a reasonable default
            for i=1:size(Rmout,2)
                Rmout[i,i]= 1.0
            end            
        end
        return Rmout
    end
    return MaterialOrientation(g!,Rm);
end



# Update the material orientation matrix.
# After the return it can be retrieved from the buffer (self.Rm).
function updateRm!(self::MaterialOrientation,XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
    self.Rm=self.g!(self.Rm, XYZ,tangents,fe_label)
    return self
end
export updateRm!

#  

function genisoRm(XYZ::JFFltMat,tangents::JFFltMat,fe_label::JFInt)
    #  Material orientation matrix for isotropic materials.
    # 
    #  function Rm = geniso_Rm(XYZ,tangents,fe_label)
    # 
    #  XYZ = location at which the material directions are needed
    #  tangents = tangent vectors to parametric coordinates in columns
    #  fe_label= label of the finite element
    # 
    #  The basic assumption here is that the material is isotropic, and
    #  therefore the choice of the material directions does not really matter as
    #  long as they correspond to the dimensionality of the element. For
    #  instance a one-dimensional element (L2 as an example) may be embedded
    #  in a three-dimensional space.
    # 
    #  This function assumes that it is being called for
    #  an ntan-dimensional manifold element, which is embedded in a
    #  sdim-dimensional Euclidean space. If ntan == sdim,
    #  the material directions matrix is the identity; otherwise the local
    #  material directions are aligned with the linear subspace defined by the
    #  tangent vectors.
    # 
    #  Warning: this *cannot* be reliably used to produce consistent stresses
    #  because each quadrature point gets a local coordinate system which
    #  depends on the orientation of the element.
    # 
    sdim, ntan = size(tangents);
    if sdim==ntan
        Rm=eye(sdim);
    else
        e1=tangents[:,1]/norm(tangents[:,1]);
        if ntan==1
            Rm = [e1];
        elseif ntan==2
            n =skewmat(e1)*tangents(:,2)/norm(tangents(:,2));
            e2=skewmat(n)*e1;
            e2=e2/norm(e2);
            Rm = [e1,e2];
        else
            error("Got an incorrect size of tangents");
        end
        return Rm
    end
end



end
