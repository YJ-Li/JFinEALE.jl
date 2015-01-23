module MaterialDeformationModule


using JFinEALE.JFFoundationModule
using JFinEALE.DeformationModelReductionModule

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


################################################################################
# 3-D model

function  stressvectorrotation{MR<:DeformationModelReduction3D}(::Type{MR},
                                                                   Rm::JFFltMat)
    # Calculate the rotation matrix for a stress vector.
    #
    #   T = stress_vector_rotation(self,Rm)
    #
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' stress vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    # Example:
    # The stress vector "stress" is given in the material coordinate
    # system defined by the orientation matrix Rm. The following two
    # transformations are equivalent:
    #
    #         t = stress_6v_to_3x3t (mat,stress);
    #         t = (Rm*t*Rm');# in global coordinate system
    #         t = (outputRm'*t*outputRm);# in output coordinate system
    #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
    #
    #        stress =mat.stress_vector_rotation(outputRm)...
    #                  *mat.stress_vector_rotation(Rm')...
    #                      *stress;# in output coordinate system
    
    
    # # Derivation of the transformation matrix [T]
    # #This is from Barbero''s  book Finite element analysis of composite
    # #materials  using Abaqus.  Note that his matrix "a"  is the transpose of
    # #the FinEALE matrix "Rm".
    # #  We also use the FinEALE numbering of the strains.
    #
    # syms T alpha R
    # syms a a11 a12 a13 a21 a22 a23 a31 a32 a33 real
    # a = [a11,a12,a13;
    #     a21,a22,a23;
    #     a31,a32,a33];
    # a = a';#his matrix "a"  is the transpose of the FinEALE matrix "Rm".
    # # # it can be done in terms of l,m,n's as well
    # # syms a l1 m1 n1 l2 m2 n2 l3 m3 n3
    # # a = [l1,m1,n1;l2,m2,n2;l3,m3,n3]
    # #  We also use the FinEALE numbering of the strains.
    # Numbering =[1,4,5;
    #     4,2,6;
    #     5,6,3];
    # T(1:6,1:6) = 0;
    # for i=1:1:3
    #     for j=1:1:3
    #         #if i==j; alpha = j; else alpha = 9-i-j; end
    #         alpha = Numbering(i,j);
    #         for p=1:1:3
    #             for q=1:1:3
    #                 #   if p==q beta = p; else beta = 9-p-q; end
    #                 beta = Numbering(p,q);
    #                 T(alpha,beta) = 0;
    #                 if alpha<=3 & beta<= 3; T(alpha,beta)=a(i,p)*a(i,p); end
    #                 if alpha> 3 & beta<= 3; T(alpha,beta)=a(i,p)*a(j,p); end
    #                 if alpha<=3 & beta>3; T(alpha,beta)=a(i,q)*a(i,p)+a(i,p)*a(i,q);end
    #                 if alpha>3 & beta>3; T(alpha,beta)=a(i,p)*a(j,q)+a(i,q)*a(j,p);end
    #             end
    #         end
    #     end
    # end
    # T
    # R = eye(6,6); R(4,4)=2; R(5,5)=2; R(6,6)=2; # Reuter matrix
    # Tbar = R*T*R^(-1)
    a11=Rm[1,1];a12=Rm[1,2];a13=Rm[1,3];
    a21=Rm[2,1];a22=Rm[2,2];a23=Rm[2,3];
    a31=Rm[3,1];a32=Rm[3,2];a33=Rm[3,3];
    T =[
        [   a11^2   a21^2   a31^2         2*a11*a21         2*a11*a31         2*a21*a31]
        [   a12^2   a22^2   a32^2         2*a12*a22         2*a12*a32         2*a22*a32]
        [   a13^2   a23^2   a33^2         2*a13*a23         2*a13*a33         2*a23*a33]
        [ a11*a12 a21*a22 a31*a32 a11*a22 + a12*a21 a11*a32 + a12*a31 a21*a32 + a22*a31]
        [ a11*a13 a21*a23 a31*a33 a11*a23 + a13*a21 a11*a33 + a13*a31 a21*a33 + a23*a31]
        [ a12*a13 a22*a23 a32*a33 a12*a23 + a13*a22 a12*a33 + a13*a32 a22*a33 + a23*a32]];
    return T
end
export stressvectorrotation

function strainvectorrotation{MR<:DeformationModelReduction3D}(::Type{MR},
                                                                      Rm::JFFltMat)
    # Calculate the rotation matrix for a strain vector.
    #
    #   Tbar = strain_vector_rotation(self,Rm)
    #
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' strain vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    
    a11=Rm[1,1];a12=Rm[1,2];a13=Rm[1,3];
    a21=Rm[2,1];a22=Rm[2,2];a23=Rm[2,3];
    a31=Rm[3,1];a32=Rm[3,2];a33=Rm[3,3];
    Tbar =[
           [     a11^2     a21^2     a31^2           a11*a21           a11*a31           a21*a31]
           [     a12^2     a22^2     a32^2           a12*a22           a12*a32           a22*a32]
           [     a13^2     a23^2     a33^2           a13*a23           a13*a33           a23*a33]
           [ 2*a11*a12 2*a21*a22 2*a31*a32 a11*a22 + a12*a21 a11*a32 + a12*a31 a21*a32 + a22*a31]
           [ 2*a11*a13 2*a21*a23 2*a31*a33 a11*a23 + a13*a21 a11*a33 + a13*a31 a21*a33 + a23*a31]
           [ 2*a12*a13 2*a22*a23 2*a32*a33 a12*a23 + a13*a22 a12*a33 + a13*a32 a22*a33 + a23*a32]];
    return Tbar
end
export strainvectorrotation

################################################################################
# 2-D plane stress model

function stressvectorrotation{MR<:DeformationModelReduction2DStress}(::Type{MR},
                                                                         Rm::JFFltMat)
    
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' stress vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    # Example:
    # The stress vector "stress" is given in the material coordinate
    # system defined by the orientation matrix Rm. The following two
    # transformations are equivalent:
    #
    #         t = stress_6v_to_3x3t (mat,stress);
    #         t = (Rm*t*Rm');# in global coordinate system
    #         t = (outputRm'*t*outputRm);# in output coordinate system
    #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
    #
    #        stress =mat.stress_vector_rotation(outputRm)...
    #                  *mat.stress_vector_rotation(Rm')...
    #                      *stress;# in output coordinate system
    
    a11=Rm[1,1];a12=Rm[1,2];
    a21=Rm[2,1];a22=Rm[2,2];
    # switch self.reduction
    #     case {'axisymm','strain'}
    #         T =[[   a11^2,   a21^2, 0,         2*a11*a21]
    #             [   a12^2,   a22^2, 0,         2*a12*a22]
    #             [       0,       0, 1,                 0]
    #             [ a11*a12, a21*a22, 0, a11*a22 + a12*a21]];
    #             case 'stress'
    T =[[   a11^2   a21^2       2*a11*a21]
        [   a12^2   a22^2       2*a12*a22]
        [ a11*a12 a21*a22 a11*a22+a12*a21]];
    return T
end
export stressvectorrotation

function strainvectorrotation{MR<:DeformationModelReduction2DStress}(::Type{MR},
                                                                            Rm::JFFltMat)
    # Calculate the rotation matrix for a strain vector.
    #
    #   Tbar = strain_vector_rotation(self,Rm)
    #
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' strain vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    
    a11=Rm[1,1];a12=Rm[1,2];
    a21=Rm[2,1];a22=Rm[2,2];
    # switch self.reduction
    #     case {'axisymm','strain'}
    #         Tbar =[ [     a11^2,     a21^2, 0,           a11*a21]
    #                 [     a12^2,     a22^2, 0,           a12*a22]
    #                 [         0,         0, 1,                 0]
    #                 [ 2*a11*a12, 2*a21*a22, 0, a11*a22 + a12*a21]];
    #   case 'stress'
    Tbar =[ [     a11^2,     a21^2,           a11*a21]
           [     a12^2,     a22^2,           a12*a22]
           [ 2*a11*a12, 2*a21*a22, a11*a22 + a12*a21]];
    return Tbar
end
export strainvectorrotation


################################################################################
# 2-D plane axially symmetric model

function stressvectorrotation{MR<:DeformationModelReduction2DAxisymm}(::Type{MR},
                                                                         Rm::JFFltMat)
    
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' stress vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    # Example:
    # The stress vector "stress" is given in the material coordinate
    # system defined by the orientation matrix Rm. The following two
    # transformations are equivalent:
    #
    #         t = stress_6v_to_3x3t (mat,stress);
    #         t = (Rm*t*Rm');# in global coordinate system
    #         t = (outputRm'*t*outputRm);# in output coordinate system
    #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
    #
    #        stress =mat.stress_vector_rotation(outputRm)...
    #                  *mat.stress_vector_rotation(Rm')...
    #                      *stress;# in output coordinate system
    
    a11=Rm[1,1];a12=Rm[1,2];
    a21=Rm[2,1];a22=Rm[2,2];
    T =[[   a11^2    a21^2  0          2*a11*a21]
        [   a12^2    a22^2  0          2*a12*a22]
        [       0        0  1                  0]
        [ a11*a12  a21*a22  0  a11*a22 + a12*a21]];
    return T
end
export stressvectorrotation

function strainvectorrotation{MR<:DeformationModelReduction2DAxisymm}(::Type{MR},
                                                                            Rm::JFFltMat)
    # Calculate the rotation matrix for a strain vector.
    #
    #   Tbar = strain_vector_rotation(self,Rm)
    #
    # Rm = columns are components of 'bar' basis vectors on the 'plain'
    #      basis vectors
    #
    # Calculate the rotation of the 'plain' strain vector to the
    # 'bar' coordinate system given by the columns of the rotation matrix Rm.
    #
    
    a11=Rm[1,1];a12=Rm[1,2];
    a21=Rm[2,1];a22=Rm[2,2];
    Tbar =[ [     a11^2      a21^2  0            a11*a21]
           [     a12^2      a22^2  0            a12*a22]
           [         0          0  1                  0]
           [ 2*a11*a12  2*a21*a22  0  a11*a22 + a12*a21]];
    return Tbar
end
export strainvectorrotation


################################################################################
# Generic versions of rotations of stiffness and compliance matrices

function rotatestiffness!{MR<:DeformationModelReduction}(::Type{MR},D::JFFltMat,
                                                           Rm::JFFltMat)
    # Rotate constitutive stiffness matrix of the material.
    #
    #         function D=transform_stiffness(self,D,Rm)
    #
    # Rotate constitutive stiffness matrix of the material to the
    # coordinate system given by the columns of the rotation matrix Rm.
    T =stressvectorrotation(MR,Rm);
    D = T*D*T';
    return D
end
export rotatestiffness!

function rotatecompliance!{MR<:DeformationModelReduction}(::Type{MR},C::JFFltMat,
                                                          Rm::JFFltMat)
    # Rotate constitutive compliance matrix of the material.
    #
    #   C = rotate_compliance(self,C,Rm)
    #
    # Rotate constitutive compliance matrix of the material to the
    # coordinate system given by the columns of the rotation matrix Rm.
    Tbar =strainvectorrotation(MR,Rm);
    C = Tbar*C*Tbar';
    return C
end
export rotatecompliance!


end
