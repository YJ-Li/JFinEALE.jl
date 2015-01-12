module FESetModule

import Base.count
import Base.cat

using JFFoundationModule
using RotationUtilModule

abstract FESet
abstract FESet0Manifold <: FESet
abstract FESet1Manifold <: FESet
abstract FESet2Manifold <: FESet
abstract FESet3Manifold <: FESet
export FESet, FESet0Manifold, FESet1Manifold, FESet2Manifold, FESet3Manifold

# Get the manifold dimension.
manifdim(me::FESet0Manifold)  =0
manifdim(me::FESet1Manifold)  =1
manifdim(me::FESet2Manifold)  =2
manifdim(me::FESet3Manifold)  =3
export manifdim

#############################################################################
# FESet
#############################################################################

# Get the number of connected finite elements
nfense{T<:FESet}(me::T) = size(me.conn,2)
export nfense

# Get the number of individual connectivities in the FE set
count{T<:FESet}(me::T) = size(me.conn,1)
export count

function getconn!{T<:FESet}(self::T,conn::JFIntMat,j::JFInt)
    # % Retrieve the connectivity for element j.
    #  
    # % Inputs//output:
    # % conn= buffer to hold the connectivity
    # %
    for i=1:size(self.conn,2)
        conn[i]=self.conn[j,i];
    end
    return self
end
export getconn!

# Evaluate the Jacobian matrix.
function Jacobianmatrix!{T<:FESet}(me::T, J::JFFltMat, x::JFFltMat, gradNparams::JFFltMat)
    for nx=1:size(x,2)
        for mx=1:size(gradNparams,2)
            J[mx,nx]=0.0;
            for kx=1:size(x,1)
                J[mx,nx]=J[mx,nx] + x[kx,mx]*gradNparams[kx,nx]; # = x' * gradNparams
            end
        end
    end
    return me
end


function setotherdimension!{T<:FESet}(self::T,other_dimension =1.0)
    if typeof(other_dimension)<:FloatingPoint
        function g(conn::JFIntMat, N::JFFltMat, x::JFFltMat)
            o::JFFlt=other_dimension::JFFlt;
            return o
        end
    else
        function g(conn::JFIntMat, N::JFFltMat, x::JFFltMat)
            o::JFFlt=other_dimension(conn,N,x)::JFFlt;
            return o
        end
    end
    self.getotherdimension=g;
    #show( code_typed(self.getotherdimension, (JFIntMat, JFFltMat, JFFltMat)) )
    return self
end
export setotherdimension!

function setlabel!{T<:FESet}(self::T,val)
    #    Set the label of this set. 
    if (length(val)==size(self.conn,1))
        self.label=val;
    elseif (length(val)==1)
        self.label=val+zeros(typeof( self.label),size(self.conn,1));
    else
        # No label intended;
        self.label=[]
    end
end
export setlabel!

function subset{T<:FESet}(self::T,L::JFIntVec)
    # Extract a subset of the individual connectivities in the FE set.
    result =deepcopy(self)
    result.conn =self.conn[L,:]
    if self.label!= nothing
        if length(self.label)==size(result.conn,1)
            result.label=self.label[L]
        end 
    end
    return  result
end
export subset

function cat{T<:FESet}(self::T, other::T)
    # % Concatenate the connectivities of two FE sets
    # %
    # % function self = cat(self, other)
    # % where
    # %  self, other = members of the same class, descendent of the FESet
    result =deepcopy(self)
    result.conn =[self.conn;other.conn];
    if (self.label!= nothing) || (other.label!= nothing)
        if (self.label== nothing)
            self.label=zeros(JFIntVec,size(self.conn,1));
        end
        if (other.label== nothing)
            other.label=zeros(JFIntVec,size(other.conn,1));
        end
    end
    result.label=[self.label;other.label];
    return result
end
export cat


function updateconn!{T<:FESet}(self::T,NewIDs::JFIntVec)
    # % Update the connectivity after the IDs of nodes changed.
    #  
    # % Inputs:
    # % NewIDs= new node IDs. Note that indexes in the conn array point
    # % _into_ the  NewIDs array. After the connectivity was updated
    # % this is no longer true!
    # %
    for i=1:size(self.conn,1)
        for j=1:size(self.conn,2)
            self.conn[i,j]=NewIDs[self.conn[i,j]];
        end        
    end
    return self
end
export updateconn!


#############################################################################
# FESet1Manifold
#############################################################################


function Jacobiancurve(self::FESet1Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the curve Jacobian.
    # 
    # %   where
    # %      self=FE set,
    # %      conn=connectivity of a single FE in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac=1.0;
    sdim, ntan = size(J);
    if     ntan==1 # 1-D FE set
        Jac = norm(J);
    else
        error("Got an incorrect size of tangents");
    end
    return Jac::JFFlt
end

function Jacobiansurface(self::FESet1Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the surface Jacobian.
    # %
    # %   For the one-dimensional cell, the surface Jacobian is 
    # %       (i) the product of the curve Jacobian and the other dimension
    # %       (units of length);
    # %       or, when used as axially symmetric
    # %       (ii) the product of the curve Jacobian and the circumference of
    # %       the circle through the point pc.
    # %    
    # %   where
    # %      self=FE set,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac=Jacobiancurve(self, conn, N, J, x);
    if self.axisymm
        xyz =N'*x;
        Jac = Jac*2*pi*xyz(1);
    else
        Jac = Jac*self.getotherdimension(conn, N, x)::JFFlt;
    end
    return Jac::JFFlt
end
        

function Jacobianvolume(self::FESet1Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the volume Jacobian.
    # %
    # %   For the one-dimensional cell, the volume Jacobian is 
    # %       (i) the product of the curve Jacobian and the other dimension
    # %       (units of length squared);
    # %       or, when used as axially symmetric
    # %       (ii) the product of the curve Jacobian and the circumference of
    # %       the circle through the point pc and the other dimension (units of
    #                                                                  %       length)
    # %    
    # %   where
    # %      self=FE set,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac= Jacobiancurve(self, conn, N, J, x);
    if self.axisymm
        xyz =N'*x;
        Jac =Jac*2*pi*xyz(1)*self.getotherdimension(conn, N, x)::JFFlt;
    else
        Jac =Jac*self.getotherdimension(conn, N, x)::JFFlt;
    end
    return Jac::JFFlt
end

function Jacobianmdim(self::FESet1Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat, m::JFInt)
    # % Evaluate the manifold Jacobian.
    # %
    # %   For and m-dimensional cell, the manifold Jacobian is
    # %       m=0: +1
    # %       m=1: Jacobian_curve
    # %       m=2: Jacobian_surface
    # %       m=3: Jacobian_volume
    # %
    # %   where
    # %      self=FE set,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac=1.0;
    if (m==3)
        Jac = Jacobianvolume(self, conn, N, J, x);
    elseif (m==2)
        Jac = Jacobiansurface(self, conn, N, J, x);
    elseif (m==1)
        Jac = Jacobiancurve(self, conn, N, J, x);
    else
        error("Wrong dimension");
    end
    return Jac::JFFlt
end

function Jacobian(self::FESet1Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the manifold Jacobian.
    # %
    # %   For and m-dimensional cell, the manifold Jacobian is
    # %       m=0: +1
    # %       m=1: Jacobian_curve
    # %       m=2: Jacobian_surface
    # %       m=3: Jacobian_volume
    # %
    # %   where
    # %      self=FE set,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac = Jacobian_curve(self, conn, N, J, x)::JFFlt;
end
        
function gradN!(self::FESet1Manifold,gradN::JFFltMat,gradNparams::JFFltMat,redJ::JFFltMat)
    # Compute the gradient of the basis functions with the respect to the "reduced" spatial coordinates.
    # 
    # gradN= output, matrix of gradients, one per row
    # gradNparams= matrix of gradients with respect to parametric coordinates, one per row
    # redJ= reduced Jacobian matrix redJ=Rm'*J
    # 
    
    for r=1:size(gradN,1)
        gradN[r,1]=  gradNparams[r,1]/redJ[1,1]
    end    
end


#############################################################################
# FESet2Manifold
#############################################################################

function Jacobiansurface(self::FESet2Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # Evaluate the surface Jacobian.
    
    #        conn=connectivity of a single cell in which the Jacobian is to be evaluated
    #        N=values of the basis functions, size(N) = [nbfuns,1]
    #        J = Jacobian matrix
    #        x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]

    Jac=1.0;
    sdim, ntan = size(J);
    if     ntan==2 # 2-D Set
        if sdim==ntan
            Jac = J[1,1]*J[2,2]-J[2,1]*J[1,2];# is det(J);% Compute the Jacobian
        else
            Jac = norm(RotationUtilModule.cross(J[:,1],J[:,2]));
        end
    else
        error("Got an incorrect size of tangents");
    end
    return Jac::JFFlt
end


function Jacobianvolume(self::FESet2Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)

    #             %   For the two-dimensional cell, the volume Jacobian is
    #             %       (i) the product of the 2-D Jacobian and the other dimension (Units
    #             %       of length);
    #             %       or, when used as axially symmetric
    #             %       (ii) the product of the 2-D Jacobian and the circumference of
    #             %       the circle through the point pc
    #             %
    #             %   where
    #             %      self=FE set,
    #             %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    #             %      N=values of the basis functions, size(N) = [nbfuns,1]
    #             %      J = Jacobian matrix
    #             %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    #             %

    Jac = Jacobiansurface(self, conn, N, J, x)::JFFlt;
    if self.axisymm
        xyz =N'*x;
        Jac = Jac*2*pi*xyz[1];
    else
        Jac = Jac*self.getotherdimension(conn, N, x)::JFFlt;
    end
    return Jac::JFFlt
end

function Jacobianmdim(self::FESet2Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat, m::JFInt)
    #             % Evaluate the manifold Jacobian.
    #             %
    #             % function Jac = Jacobian_mdim(self, conn, N, J, x, m)
    #             %
    #             %
    #             %   For and m-dimensional cell, the manifold Jacobian is
    #             %       m=0: +1
    #             %       m=1: Jacobian_curve
    #             %       m=2: Jacobian_surface
    #             %       m=3: Jacobian_volume
    #             %
    #             %   where
    #             %      self=FE set,
    #             %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    #             %      N=values of the basis functions, size(N) = [nbfuns,1]
    #             %      J = Jacobian matrix
    #             %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    #             %
    Jac=1.0;
    if m==3
        Jac = Jacobianvolume(self, conn, N, J, x)::JFFlt;
    elseif m==2
        Jac = Jacobiansurface(self, conn, N, J, x)::JFFlt;
    else
        error("Wrong dimension"); Jac::JFFlt=0.0;
    end
    return Jac::JFFlt
end

function Jacobian(self::FESet2Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    #   % Evaluate the manifold Jacobian.
    #   %
    #   % function Jac = Jacobian(self, conn, N, J, x)
    #     %
    #     %   For and m-dimensional cell, the manifold Jacobian is
    #     %       m=0: +1
    #     %       m=1: Jacobian_curve
    #     %       m=2: Jacobian_surface
    #     %       m=3: Jacobian_volume
    #     %
    #     %   where
    #     %      self=FE set,
    #     %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    #     %      N=values of the basis functions, size(N) = [nbfuns,1]
    #     %      J = Jacobian matrix
    #     %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    #     %
    Jac = Jacobiansurface(self, conn, N, J, x)::JFFlt;
end

function gradN!(self::FESet2Manifold,gradN::JFFltMat,gradNparams::JFFltMat,redJ::JFFltMat)
    # Compute the gradient of the basis functions with the respect to the "reduced" spatial coordinates.
    # 
    # gradN= output, matrix of gradients, one per row
    # gradNparams= matrix of gradients with respect to parametric coordinates, one per row
    # redJ= reduced Jacobian matrix redJ=Rm'*J
    #
    
    # detredJ=(redJ[1,1]*redJ[2,2] - redJ[1,2]*redJ[2,1])
    # ext=[redJ[2,2]/detredJ -redJ[1,2]/detredJ;
    #      -redJ[2,1]/detredJ  redJ[1,1]/detredJ]
    # for c=1:size(gradN,2)
    #     for r=1:size(gradN,1)
    #         gradN[r,c]=0.0
    #         for k=1:size(gradN,2)
    #             gradN[r,c]= gradN[r,c]+gradNparams[r,k]*ext[k,c]
    #         end
    #     end
    # end
    
    # This is the unrolled version that avoids allocation of a 3 x 3 matrix
    invdet=1.0/(redJ[1,1]*redJ[2,2] - redJ[1,2]*redJ[2,1]);
    invredJ11 =  (redJ[2,2])*invdet;
    invredJ12 = -(redJ[1,2])*invdet;
    invredJ21 = -(redJ[2,1])*invdet;
    invredJ22 =  (redJ[1,1])*invdet;
    for r=1:size(gradN,1)
        gradN[r,1]= gradNparams[r,1]*invredJ11 +gradNparams[r,2]*invredJ21;
        gradN[r,2]= gradNparams[r,1]*invredJ12 +gradNparams[r,2]*invredJ22;
    end
end


#############################################################################
# FESet3Manifold
#############################################################################


function Jacobianvolume(self::FESet3Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the volume Jacobian.
    # %
    # %   where
    # %      self=gcellset,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac=1.0;
    sdim, ntan = size(J);
    if     ntan==3
        Jac = det(J);# Compute the Jacobian
    else
        error("Got an incorrect size of tangents");
    end
    return Jac::JFFlt
end

function Jacobianmdim(self::FESet3Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat, m::JFInt)
    # % Evaluate the manifold Jacobian.
    # %
    # %   For and m-dimensional cell, the manifold Jacobian is
    # %       m=0: +1
    # %       m=1: Jacobian_curve
    # %       m=2: Jacobian_surface
    # %       m=3: Jacobian_volume
    # %
    # %   where
    # %      self=gcellset,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac=1.0;
    if (m==3)
        Jac = Jacobianvolume(self, conn, N, J, x);
    else
        error("Wrong dimension"); 
    end
    return Jac::JFFlt
end

function Jacobian(self::FESet3Manifold, conn::JFIntMat, N::JFFltMat, J::JFFltMat, x::JFFltMat)
    # % Evaluate the manifold Jacobian.
    # %
    # %   For and m-dimensional cell, the manifold Jacobian is
    # %       m=0: +1
    # %       m=1: Jacobian_curve
    # %       m=2: Jacobian_surface
    # %       m=3: Jacobian_volume
    # %
    # %   where
    # %      self=gcellset,
    # %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
    # %      N=values of the basis functions, size(N) = [nbfuns,1]
    # %      J = Jacobian matrix
    # %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
    # %
    Jac = Jacobian_volume(self, conn, N, J, x)::JFFlt;
end
        
function gradN!(self::FESet3Manifold,gradN::JFFltMat,gradNparams::JFFltMat,redJ::JFFltMat)
    # Compute the gradient of the basis functions with the respect to the "reduced" spatial coordinates.
    # 
    # gradN= output, matrix of gradients, one per row
    # gradNparams= matrix of gradients with respect to parametric coordinates, one per row
    # redJ= reduced Jacobian matrix redJ=Rm'*J
    # 
    invdet = 1.0 / ( +redJ[1,1]*(redJ[2,2]*redJ[3,3]-redJ[3,2]*redJ[2,3])
                    -redJ[1,2]*(redJ[2,1]*redJ[3,3]-redJ[2,3]*redJ[3,1])
                    +redJ[1,3]*(redJ[2,1]*redJ[3,2]-redJ[2,2]*redJ[3,1]) );
    # local invredJ=zeros(JFFlt,3,3);
    # invredJ[1,1] =  (redJ[2,2]*redJ[3,3]-redJ[3,2]*redJ[2,3])*invdet;
    # invredJ[1,2] = -(redJ[1,2]*redJ[3,3]-redJ[1,3]*redJ[3,2])*invdet;
    # invredJ[1,3] =  (redJ[1,2]*redJ[2,3]-redJ[1,3]*redJ[2,2])*invdet;
    # invredJ[2,1] = -(redJ[2,1]*redJ[3,3]-redJ[2,3]*redJ[3,1])*invdet;
    # invredJ[2,2] =  (redJ[1,1]*redJ[3,3]-redJ[1,3]*redJ[3,1])*invdet;
    # invredJ[2,3] = -(redJ[1,1]*redJ[2,3]-redJ[2,1]*redJ[1,3])*invdet;
    # invredJ[3,1] =  (redJ[2,1]*redJ[3,2]-redJ[3,1]*redJ[2,2])*invdet;
    # invredJ[3,2] = -(redJ[1,1]*redJ[3,2]-redJ[3,1]*redJ[1,2])*invdet;
    # invredJ[3,3] =  (redJ[1,1]*redJ[2,2]-redJ[2,1]*redJ[1,2])*invdet;
    # for c=1:size(gradN,2)
    #     for r=1:size(gradN,1)
    #         gradN[r,c]=0.0
    #         for k=1:size(gradN,2)
    #             gradN[r,c]= gradN[r,c]+gradNparams[r,k]*invredJ[k,c]
    #         end
    #     end
    # end

    # This is the unrolled version that avoids allocation of a 3 x 3 matrix
    invredJ11 =  (redJ[2,2]*redJ[3,3]-redJ[3,2]*redJ[2,3])*invdet;
    invredJ12 = -(redJ[1,2]*redJ[3,3]-redJ[1,3]*redJ[3,2])*invdet;
    invredJ13 =  (redJ[1,2]*redJ[2,3]-redJ[1,3]*redJ[2,2])*invdet;
    invredJ21 = -(redJ[2,1]*redJ[3,3]-redJ[2,3]*redJ[3,1])*invdet;
    invredJ22 =  (redJ[1,1]*redJ[3,3]-redJ[1,3]*redJ[3,1])*invdet;
    invredJ23 = -(redJ[1,1]*redJ[2,3]-redJ[2,1]*redJ[1,3])*invdet;
    invredJ31 =  (redJ[2,1]*redJ[3,2]-redJ[3,1]*redJ[2,2])*invdet;
    invredJ32 = -(redJ[1,1]*redJ[3,2]-redJ[3,1]*redJ[1,2])*invdet;
    invredJ33 =  (redJ[1,1]*redJ[2,2]-redJ[2,1]*redJ[1,2])*invdet;
    for r=1:size(gradN,1)
        gradN[r,1]= gradNparams[r,1]*invredJ11 +gradNparams[r,2]*invredJ21 +gradNparams[r,3]*invredJ31;
        gradN[r,2]= gradNparams[r,1]*invredJ12 +gradNparams[r,2]*invredJ22 +gradNparams[r,3]*invredJ32;
        gradN[r,3]= gradNparams[r,1]*invredJ13 +gradNparams[r,2]*invredJ23 +gradNparams[r,3]*invredJ33;
    end
end


#############################################################################
# FESetP1
#############################################################################

type FESetP1 <: FESet0Manifold
    #        # Connectivity array of the set.
    #        #    The connectivity array lists the nodes that the finite element in
    #        #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetP1(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=1)
            error("Wrong number of nodes.");
        end
        # Need to make a COPY of the input arrays
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end
export FESetP1

bfun(self::FESetP1, param_coords::JFFltMat) = reshape([1.0],1,1) # make sure this is a matrix
bfundpar(self::FESetP1, param_coords::JFFltMat) = reshape([0.0],1,1)

function boundaryconn(self::FESetP1)
    # Get boundary connectivity.
    return [];
end

function boundaryfe(self::FESetP1)
    # Get  the constructor of the class of the  boundary finite element.
    return None;
end

#############################################################################
# FESetL2
#############################################################################

type FESetL2 <: FESet1Manifold
    #        # Connectivity array of the set.
    #        #    The connectivity array lists the nodes that the finite element in
    #        #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetL2(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=2)
            error("Wrong number of nodes.");
        end
        # Needs to make a COPY of the input arrays
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end
export FESetL2 

bfun(self::FESetL2, param_coords::JFFltMat) = reshape([(1. - param_coords[1]); (1. + param_coords[1])] / 2.0,2,1) # make sure this is a matrix
bfundpar(self::FESetL2, param_coords::JFFltMat) = reshape([-1.0; +1.0]/2.0,2,1)

function boundaryconn(self::FESetL2)
    # Get boundary connectivity.
    return [self.conn[:,1];self.conn[:,2]];
end

function boundaryfe(self::FESetL2)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetP1;
end

#############################################################################
# FESetT3
#############################################################################

type FESetT3 <: FESet2Manifold
    # Connectivity array of the set.
    #    The connectivity array lists the nodes that the finite element in
    #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetT3(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=3)
            error("Wrong number of nodes.");
        end
        # Needs to make a COPY of the input arrays
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end
export FESetT3

function bfun(self::FESetT3, param_coords::JFFltMat)
    # Evaluate the basis function matrix for an 3-node triangle.
    return reshape([(1 - param_coords[1] - param_coords[2]); param_coords[1]; param_coords[2]],3,1); # Make sure this is a matrix
end

function bfundpar(self::FESetT3, param_coords::JFFltMat)
    # Evaluate the derivatives of the basis function matrix.
    return  [-1. -1.;  +1.  0.;  0. +1.];
end

function boundaryconn(self::FESetT3)
    # Get boundary connectivity.
    return [self.conn[:,1:2];self.conn[:,2:3];self.conn[:,[3,1]]];
end

function boundaryfe(self::FESetT3)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end


#############################################################################
# FESetQ4
#############################################################################

type FESetQ4 <: FESet2Manifold
    # Connectivity array of the set.
    #    The connectivity array lists the nodes that the finite element in
    #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetQ4(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=4)
            error("Wrong number of nodes.");
        end
        # Needs to make a COPY of the input arrays
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end
export FESetQ4 

function bfun(self::FESetQ4, param_coords::JFFltMat)
    # Evaluate the basis function matrix for an 4-node quadrilateral.
    val = [0.25 * (1. - param_coords[1]) * (1. - param_coords[2]);
           0.25 * (1. + param_coords[1]) * (1. - param_coords[2]);
           0.25 * (1. + param_coords[1]) * (1. + param_coords[2]);
           0.25 * (1. - param_coords[1]) * (1. + param_coords[2])];
    return reshape(val,4,1); # Make sure this is a matrix
end

function bfundpar(self::FESetQ4, param_coords::JFFltMat)
    # Evaluate the derivatives of the basis function matrix.
    val = [-(1. - param_coords[2])*0.25 -(1. - param_coords[1])*0.25;
            (1. - param_coords[2])*0.25 -(1. + param_coords[1])*0.25;
            (1. + param_coords[2])*0.25 (1. + param_coords[1])*0.25;
           -(1. + param_coords[2])*0.25 (1. - param_coords[1])*0.25];
    return val
end

function boundaryconn(self::FESetQ4)
    # Get boundary connectivity.
    return [self.conn[:,1:2];self.conn[:,2:3];self.conn[:,3:4];self.conn[:,[4,1]]];
end

function boundaryfe(self::FESetQ4)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetL2;
end


#############################################################################
# FESetH8
#############################################################################

type FESetH8 <: FESet3Manifold
    # Connectivity array of the set.
    #    The connectivity array lists the nodes that the finite element in
    #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetH8(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=8)
            error("Wrong number of nodes.");
        end
        # Needs to make a COPY of the input arrays
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end
export FESetH8

function bfun(self::FESetH8, param_coords::JFFltMat)
    # Evaluate the basis function matrix for an 8-node hexahedron.
    one_minus_xi    = (1.0 - param_coords[1]);
    one_minus_eta   = (1.0 - param_coords[2]);
    one_minus_theta = (1.0 - param_coords[3]);
    one_plus_xi     = (1.0 + param_coords[1]);
    one_plus_eta    = (1.0 + param_coords[2]);
    one_plus_theta  = (1.0 + param_coords[3]);
    val = [one_minus_xi*one_minus_eta*one_minus_theta;
           one_plus_xi*one_minus_eta*one_minus_theta;
           one_plus_xi*one_plus_eta*one_minus_theta;
           one_minus_xi*one_plus_eta*one_minus_theta;
           one_minus_xi*one_minus_eta*one_plus_theta;
           one_plus_xi*one_minus_eta*one_plus_theta;
           one_plus_xi*one_plus_eta*one_plus_theta;
           one_minus_xi*one_plus_eta*one_plus_theta] / 8.0;
    return reshape(val,8,1); # Make sure this is a matrix
end

function bfundpar(self::FESetH8, param_coords::JFFltMat)
    # Evaluate the derivatives of the basis function matrix.
     omxi    = (1.0 - param_coords[1]);
    ometa   = (1.0 - param_coords[2]);
    omtheta = (1.0 - param_coords[3]);
    opxi     = (1.0 + param_coords[1]);
    opeta    = (1.0 + param_coords[2]);
    optheta  = (1.0 + param_coords[3]);
   val = [-ometa*omtheta  ometa*omtheta  opeta*omtheta  -opeta*omtheta  -ometa*optheta  ometa*optheta  opeta*optheta   -opeta*optheta;
          -omxi*omtheta  -opxi*omtheta   opxi*omtheta   omxi*omtheta    -omxi*optheta   -opxi*optheta  opxi*optheta   omxi*optheta;
          -omxi*ometa    -opxi*ometa     -opxi*opeta    -omxi*opeta      omxi*ometa     opxi*ometa     opxi*opeta    omxi*opeta]'/8.0;
    return val
end

function boundaryconn(self::FESetH8)
    # Get boundary connectivity.
    return  [self.conn[:,[1, 4, 3, 2]];
             self.conn[:,[1, 2, 6, 5]];
             self.conn[:,[2, 3, 7, 6]];
             self.conn[:,[3, 4, 8, 7]];
             self.conn[:,[4, 1, 5, 8]];
             self.conn[:,[6, 7, 8, 5]]];
end

function boundaryfe(self::FESetH8)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ4;
end


#############################################################################
# FESetH27
#############################################################################

type FESetH27 <: FESet3Manifold
    # Connectivity array of the set.
    #    The connectivity array lists the nodes that the finite element in
    #  the set connects. The j-th element connects the nodes conn(j,:)
    conn::JFIntMat
    #        # Numerical label, supplied for each element in the set, or
    #        # a single number to be applied to all elements
    getotherdimension::Function
    axisymm::Bool
    label::JFIntVec

    function    FESetH27(;conn::JFIntMat=[],other_dimension =1.0,axisymm::Bool=false,label =[])
        if (size(conn,2)!=27)
            error("Wrong number of nodes.");
        end
        self =new(deepcopy(conn),(x,y,z)->  1.0,axisymm,deepcopy(label))
        setotherdimension!(self,other_dimension)
        setlabel!(self,label)
        return self
    end
end

function bfun(self::FESetH27, param_coords::JFFltMat)
    # Evaluate the basis function matrix for an 8-node hexahedron.
    xi=param_coords[1];
    eta=param_coords[2];
    zet=param_coords[3];
    x1 =(xi-1);
    y1=(eta-1);
    z1 =(zet-1);
    x2 =(xi+1);
    y2=(eta+1);
    z2 =(zet+1);
    val =[1.0/8.0*z1*zet*x1*xi*y1*eta
          1.0/8.0*z1*zet*x2*xi*y1*eta
          1.0/8.0*z1*zet*x2*xi*y2*eta
          1.0/8.0*z1*zet*x1*xi*y2*eta
          1.0/8.0*z2*zet*x1*xi*y1*eta
          1.0/8.0*z2*zet*x2*xi*y1*eta
          1.0/8.0*z2*zet*x2*xi*y2*eta
          1.0/8.0*z2*zet*x1*xi*y2*eta
          1.0/4.0*z1*zet*(-x2)*x1*y1*eta
          1.0/4.0*z1*zet*x2*xi*(-y2)*y1
          1.0/4.0*z1*zet*(-x2)*x1*y2*eta
          1.0/4.0*z1*zet*x1*xi*(-y2)*y1
          1.0/4.0*z2*zet*(-x2)*x1*y1*eta
          1.0/4.0*z2*zet*x2*xi*(-y2)*y1
          1.0/4.0*z2*zet*(-x2)*x1*y2*eta
          1.0/4.0*z2*zet*x1*xi*(-y2)*y1
          1.0/4.0*(-z2)*z1*x1*xi*y1*eta
          1.0/4.0*(-z2)*z1*x2*xi*y1*eta
          1.0/4.0*(-z2)*z1*x2*xi*y2*eta
          1.0/4.0*(-z2)*z1*x1*xi*y2*eta
          1.0/2.0*z1*zet*(-x2)*x1*(-y2)*y1
          1.0/2.0*(-z2)*z1*(-x2)*x1*y1*eta
          1.0/2.0*(-z2)*z1*x2*xi*(-y2)*y1
          1.0/2.0*(-z2)*z1*(-x2)*x1*y2*eta
          1.0/2.0*(-z2)*z1*x1*xi*(-y2)*y1
          1.0/2.0*z2*zet*(-x2)*x1*(-y2)*y1
          (-z2)*z1*(-x2)*x1*(-y2)*y1];
    return reshape(val,27,1); # Make sure this is a matrix
end

function bfundpar(self::FESetH27, param_coords::JFFltMat)
    # Evaluate the derivatives of the basis function matrix.
    xi=param_coords[1];
    eta=param_coords[2];
    zet=param_coords[3];
    x1 =(xi-1.0/2.0);
    x2 =(xi+1.0/2.0);
    x3 =(xi-1.0);
    x4 =(xi+1.0);
    z1 =(zet-1.0);
    z2 =(zet-1.0/2.0);
    z3 =(zet+1.0);
    z4 =(zet+1.0/2.0);
    y1 =(eta-1.0);
    y2 =(eta-1.0/2.0);
    y3 =(eta+1.0);
    y4 = (eta+1.0/2.0);
    val = [      1.0/4.0*z1*zet*x1*y1*eta        1.0/4.0*z1*zet*x3*xi*y2        1.0/4.0*z2*x3*xi*y1*eta;
                  1.0/4.0*z1*zet*x2*y1*eta        1.0/4.0*z1*zet*x4*xi*y2        1.0/4.0*z2*x4*xi*y1*eta;
                  1.0/4.0*z1*zet*x2*y3*eta        1.0/4.0*z1*zet*x4*xi*y4        1.0/4.0*z2*x4*xi*y3*eta;
                  1.0/4.0*z1*zet*x1*y3*eta        1.0/4.0*z1*zet*x3*xi*y4        1.0/4.0*z2*x3*xi*y3*eta;
                  1.0/4.0*z3*zet*x1*y1*eta        1.0/4.0*z3*zet*x3*xi*y2        1.0/4.0*z4*x3*xi*y1*eta;
                  1.0/4.0*z3*zet*x2*y1*eta        1.0/4.0*z3*zet*x4*xi*y2        1.0/4.0*z4*x4*xi*y1*eta;
                  1.0/4.0*z3*zet*x2*y3*eta        1.0/4.0*z3*zet*x4*xi*y4        1.0/4.0*z4*x4*xi*y3*eta;
                  1.0/4.0*z3*zet*x1*y3*eta        1.0/4.0*z3*zet*x3*xi*y4        1.0/4.0*z4*x3*xi*y3*eta;
                       -1.0/2.0*z1*zet*xi*y1*eta   1.0/2.0*z1*zet*(-x4)*x3*y2   1.0/2.0*z2*(-x4)*x3*y1*eta;
             1.0/2.0*z1*zet*x2*(-y3)*y1             -1.0/2.0*z1*zet*x4*xi*eta   1.0/2.0*z2*x4*xi*(-y3)*y1;
                       -1.0/2.0*z1*zet*xi*y3*eta   1.0/2.0*z1*zet*(-x4)*x3*y4   1.0/2.0*z2*(-x4)*x3*y3*eta;
             1.0/2.0*z1*zet*x1*(-y3)*y1             -1.0/2.0*z1*zet*x3*xi*eta   1.0/2.0*z2*x3*xi*(-y3)*y1;
                       -1.0/2.0*z3*zet*xi*y1*eta   1.0/2.0*z3*zet*(-x4)*x3*y2   1.0/2.0*z4*(-x4)*x3*y1*eta;
             1.0/2.0*z3*zet*x2*(-y3)*y1             -1.0/2.0*z3*zet*x4*xi*eta   1.0/2.0*z4*x4*xi*(-y3)*y1;
                       -1.0/2.0*z3*zet*xi*y3*eta   1.0/2.0*z3*zet*(-x4)*x3*y4   1.0/2.0*z4*(-x4)*x3*y3*eta;
             1.0/2.0*z3*zet*x1*(-y3)*y1             -1.0/2.0*z3*zet*x3*xi*eta   1.0/2.0*z4*x3*xi*(-y3)*y1;
             1.0/2.0*(-z3)*z1*x1*y1*eta            1.0/2.0*(-z3)*z1*x3*xi*y2             -1.0/2.0*zet*x3*xi*y1*eta;
             1.0/2.0*(-z3)*z1*x2*y1*eta            1.0/2.0*(-z3)*z1*x4*xi*y2             -1.0/2.0*zet*x4*xi*y1*eta;
             1.0/2.0*(-z3)*z1*x2*y3*eta            1.0/2.0*(-z3)*z1*x4*xi*y4             -1.0/2.0*zet*x4*xi*y3*eta;
             1.0/2.0*(-z3)*z1*x1*y3*eta            1.0/2.0*(-z3)*z1*x3*xi*y4             -1.0/2.0*zet*x3*xi*y3*eta;
                      -z1*zet*xi*(-y3)*y1            -z1*zet*(-x4)*x3*eta  z2*(-x4)*x3*(-y3)*y1;
                      -(-z3)*z1*xi*y1*eta  (-z3)*z1*(-x4)*x3*y2            -zet*(-x4)*x3*y1*eta;
            (-z3)*z1*x2*(-y3)*y1            -(-z3)*z1*x4*xi*eta            -zet*x4*xi*(-y3)*y1;
                      -(-z3)*z1*xi*y3*eta  (-z3)*z1*(-x4)*x3*y4            -zet*(-x4)*x3*y3*eta;
            (-z3)*z1*x1*(-y3)*y1            -(-z3)*z1*x3*xi*eta            -zet*x3*xi*(-y3)*y1;
                      -z3*zet*xi*(-y3)*y1            -z3*zet*(-x4)*x3*eta  z4*(-x4)*x3*(-y3)*y1;
               -2.0*(-z3)*z1*xi*(-y3)*y1     -2.0*(-z3)*z1*(-x4)*x3*eta     -2.0*zet*(-x4)*x3*(-y3)*y1;
           ];
    return reshape(val,27,3);
end

function boundaryconn(self::FESetH27)
    # Get boundary connectivity.
    return  [self.conn[:,[1, 4, 3, 2, 12, 11, 10, 9, 21]];
             self.conn[:,[1, 2, 6, 5, 9, 18, 13, 17, 22]];
             self.conn[:,[2, 3, 7, 6, 10, 19, 14, 18, 23]];
             self.conn[:,[3, 4, 8, 7,  11, 20, 15, 19, 24]];
             self.conn[:,[4, 1, 5, 8, 12, 17, 16, 20, 25]];
             self.conn[:,[6, 7, 8, 5, 14, 15, 16, 13, 26]]];
end

function boundaryfe(self::FESetH27)
    # Get  the constructor of the class of the  boundary finite element.
    return FESetQ9;
end

end

