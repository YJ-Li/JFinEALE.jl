module MaterialDeformationModule


using JFinEALE.JFFoundationModule


# % Class that represents deformable materials.
    # %
    # % The purpose of the class is to provide a variety of conversion methods for
    # % stress and strain to all descendents of this class.


function strain2x2tto3v! (v::JFFltVec,t::JFFltMat)
    # % Convert a matrix of 2x2 strain components  into a 3-component vector.
    #  
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[1,2] + t[2,1];
end
export strain2x2tto3v!

function strain3vto2x2t! (t::JFFltMat,v::JFFltVec)
    # % Convert a strain 3-vector to a  matrix of 2x2 strain components (symmetric tensor)
    # 
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3]/2;
    t[2,1] = v[3]/2;
end
export strain3vto2x2t!

function strain3x3tto6v! (v::JFFltVec,t::JFFltMat)
    # Convert a matrix of 3x3 strain components to a 6-component vector.
    
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[3,3];
    v[4] = t[1,2] + t[2,1];
    v[5] = t[1,3] + t[3,1];
    v[6] = t[3,2] + t[2,3];
end
export strain3x3tto6v!

function strain6vto3x3t! (t::JFFltMat,v::JFFltVec)
    # Convert a strain 6-vector to a  matrix of 3x3 strain components (symmetric tensor)
    
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[3,3] = v[3];
    t[1,2] = v[4]/2.;
    t[2,1] = v[4]/2.;
    t[1,3] = v[5]/2.;
    t[3,1] = v[5]/2.;
    t[3,2] = v[6]/2.;
    t[2,3] = v[6]/2.;
end
export strain6vto3x3t!


function stress2x2to3v! (v::JFFltVec,t::JFFltMat)
    # Convert a symmetric matrix of 2x2 stress components to a 3-component vector.
    #
    
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = 1./2.*(t[1,2] + t[2,1]);
end
export stress2x2to3v!

function stress3vto2x2t! (t::JFFltMat,v::JFFltVec)
    # Convert a 3-vector to a  matrix of 2x2 stress components (symmetric tensor)
     
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
end
export stress3vto2x2t!

function stress3vto3x3t! (t::JFFltMat,v::JFFltVec)
    # Convert a 3-vector to a matrix of 3x3 stress components (symmetric tensor)
    #
    
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
end
export stress3vto3x3t!

function stress4vto3x3t! (t::JFFltMat,v::JFFltVec)
    # Convert a 4-vector to a  matrix of 3x3 stress components (tensor).
    #
    # function t = stress4vto3x3t (self,v)
    #
    #
    # Convert a 4-vector to a *symmetric*
    # matrix of 3x3 stress components (tensor).  This is
    # conversion routine that would be useful for plane strain or 
    # axially symmetric conditions.
    # The stress vector components need to be ordered as:
    #     sigmax, sigmay, tauxy, sigmaz,
    # which is the ordering used for the plane-strain model reduction.
    # Therefore, for axially symmetric analysis the components need to be
    # reordered, as from the constitutive equation they come out 
    # as sigmax, sigmay, sigmaz, tauxy.
    #
    
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    t[3,3] = v[4];
end
export stress4vto3x3t!

function stress6vto3x3t! (t::JFFltMat,v::JFFltVec)
    # Convert a 6-vector to a  matrix of 3x3 stress components (symmetric tensor)
    # 
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[3,3] = v[3];
    t[1,2] = v[4];
    t[2,1] = v[4];
    t[1,3] = v[5];
    t[3,1] = v[5];
    t[3,2] = v[6];
    t[2,3] = v[6];
end
export stress6vto3x3t!


function stress3x3tto6v! (v::JFFltVec,t::JFFltMat)
    # Convert a matrix of 3x3 stress components to a 6-component vector.
    #
    # function v = stress3x3tto6v (self,t)
    #
    # Convert a matrix of 3x3 stress components (tensor)
    # into a 6-component vector.
    #
   
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[3,3];
    v[4] = 1/2.*(t[1,2] + t[2,1]);
    v[5] = 1/2.*(t[1,3] + t[3,1]);
    v[6] = 1/2.*(t[3,2] + t[2,3]);
   # v = t[[1,5,9,2,3,6]]';
end
export stress3x3tto6v!

function strain9vto6v! (t::JFFltVec,v::JFFltVec)
    # Convert a strain 9-vector to a  strain 6-vector components (tensor)
    #
    # function t = strain9vto6v (self,v)
    #
    
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4]+v[5];
    t[5] = v[8]+v[9];
    t[6] = v[6]+v[7];
end
export strain9vto6v!

function strain6vto9v! (t::JFFltVec,v::JFFltVec)
    # Convert a strain 6-vector to a  strain 9-vector components (tensor)
    
    
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4]/2.;
    t[5] = v[4]/2.;
    t[6] = v[6]/2.;
    t[7] = v[6]/2.;
    t[8] = v[5]/2.;
    t[9] = v[5]/2.;
end
export strain6vto9v!

function stress9vto6v! (t::JFFltVec,v::JFFltVec)
    # Convert a stress 9-vector to a  stress 6-vector components (tensor)
    
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4];
    t[5] = v(8);
    t[6] = v[6];
end
export stress9vto6v!

end
