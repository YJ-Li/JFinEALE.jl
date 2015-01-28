module DeformationModelReductionModule

using JFinEALE.JFFoundationModule


abstract DeformationModelReduction;
type DeformationModelReduction1D <: DeformationModelReduction
end
type DeformationModelReduction2DStrain <: DeformationModelReduction
end
type DeformationModelReduction2DStress <: DeformationModelReduction
end
type DeformationModelReduction2DAxisymm <: DeformationModelReduction
end
type DeformationModelReduction3D <: DeformationModelReduction
end
export DeformationModelReduction
export DeformationModelReduction1D
export DeformationModelReduction2DStrain
export DeformationModelReduction2DStress
export DeformationModelReduction2DAxisymm
export DeformationModelReduction3D


nstrains(::Type{DeformationModelReduction1D})=1;
export nstrains
nstrains(::Type{DeformationModelReduction2DStrain})=3;
export nstrains
nstrains(::Type{DeformationModelReduction2DStress})=3;
export nstrains
nstrains(::Type{DeformationModelReduction2DAxisymm})=4;
export nstrains
nstrains(::Type{DeformationModelReduction3D})=6;
export nstrains

function stresscomponentmap(::Type{DeformationModelReduction1D})
    a=Dict{Symbol,Int64}(); a[:x]=1;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeformationModelReduction2DStrain})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:xy]=3;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeformationModelReduction2DStress})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:xy]=3;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeformationModelReduction2DAxisymm})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4;
    return a
end
export stresscomponentmap

function stresscomponentmap(::Type{DeformationModelReduction3D})
    a=Dict{Symbol,Int64}(); a[:x]=1;  a[:y]=2;  a[:z]=3;  a[:xy]=4; a[:xz]=5; a[:yz]=6;
    return a
end
export stresscomponentmap

function Blmat!(::Type{DeformationModelReduction1D},B::JFFltMat,N::JFFltMat,gradN::JFFltMat,c::JFFltMat,Rm::JFFltMat)
    # Compute the strain-displacement matrix for a one-manifold element.
    #
    # B = Blmat1(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a one-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [1,nfens*ndim],
    #         where  nfens = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(1,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        @inbounds for j=1:dim
            k=dim*(i-1);
            B[1,k+j]=  gradN[i,1] *Rm[j,1];
        end
    end
    return B::JFFltMat;
end
export Blmat!

function Blmat!(::Type{DeformationModelReduction2DStrain},
                B::JFFltMat,N::JFFltMat,gradN::JFFltMat,c::JFFltMat,Rm::JFFltMat)
    # Compute the strain-displacement matrix for a two-manifold element.
    #
    # B = Blmat2(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a two-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [3,nfens*2],
    #         where  nfens = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(3,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            B[1,k+j]=gradN[i,1]*Rm[j,1]
            B[2,k+j]=gradN[i,2]*Rm[j,2]
            B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
        end
    end
    return B::JFFltMat;
end
export Blmat!

function Blmat!(::Type{DeformationModelReduction2DStress},
                B::JFFltMat,N::JFFltMat,gradN::JFFltMat,c::JFFltMat,Rm::JFFltMat)
    # Compute the strain-displacement matrix for a two-manifold element.
    #
    # B = Blmat2(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a two-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [3,nfens*2],
    #         where  nfens = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(3,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            B[1,k+j]=gradN[i,1]*Rm[j,1]
            B[2,k+j]=gradN[i,2]*Rm[j,2]
            B[3,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
        end
    end
    return B::JFFltMat;
end
export Blmat!

function Blmat!(::Type{DeformationModelReduction2DAxisymm},
                B::JFFltMat,N::JFFltMat,gradN::JFFltMat,c::JFFltMat,Rm::JFFltMat)
    # Compute the strain-displacement matrix for a two-manifold element for
    # axially symmetric models.
    #
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
    #         of the embedding space, greater than or equal to  2.
    # Output:
    #      B = strain-displacement matrix, size (B) = [4,nfens*2],
    #         where  nfens = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    r=c[1];
    if r==<eps(1.0)
        r=eps(1.0);
    end
    dim =size(c,2);
    @assert size(B)==(4,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            B[1,k+j]=gradN[i,1]*Rm[j,1]
            B[2,k+j]=gradN[i,2]*Rm[j,2]
            B[3,k+j]=N[i]/r*Rm[j,1]
            B[4,k+j]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
        end
    end
   return B::JFFltMat;
end
export Blmat!

function Blmat!(::Type{DeformationModelReduction3D},
                B::JFFltMat,N::JFFltMat,gradN::JFFltMat,c::JFFltMat,Rm::JFFltMat)
    # Compute the strain-displacement matrix for a three-manifold element.
    #
    # B = Blmat3(self,N,gradN,c,Rm)
    #
    # Compute the linear, displacement independent, strain-displacement matrix
    # for a three-manifold element.
    #
    # Arguments
    #      self=finite element model
    #      N= matrix of basis function values
    #      gradN=matrix of basis function gradients with respect to the
    #         Cartesian coordinates in the directions of the material orientation
    #      c=array of spatial coordinates of the evaluation point
    #         in the global Cartesian coordinates
    #      Rm=orthogonal matrix with the unit basis vectors of the local
    #         material orientation coordinate system as columns;
    #         supply Rm as empty when there is no need for global-to-local
    #         transformation (for instance for isotropic materials).
    #         size(Rm)= [ndim,ndim], where ndim = number of spatial dimensions
    #         of the embedding space.
    # Output:
    #      B = strain-displacement matrix, size (B) = [nstrain,nfens*dim],
    #         where nstrain= number of strains, dim = Number of spatial
    #         dimensions of the embedding space, and nfens = number of
    #         finite element nodes on the element.
    #
    nfn= size(gradN,1);
    dim =size(c,2);
    @assert size(B)==(6,dim*nfn)
    fill!(B,0.0);
    @inbounds for i= 1:nfn
        k=dim*(i-1);
        @inbounds for j=1:dim
            kj=k+j
            B[1,kj]=gradN[i,1]*Rm[j,1]
            B[2,kj]=gradN[i,2]*Rm[j,2]
            B[3,kj]=gradN[i,3]*Rm[j,3]
            B[4,kj]=gradN[i,2]*Rm[j,1]+gradN[i,1]*Rm[j,2]
            B[5,kj]=gradN[i,3]*Rm[j,1]+gradN[i,1]*Rm[j,3]
            B[6,kj]=gradN[i,3]*Rm[j,2]+gradN[i,2]*Rm[j,3]
        end
    end
    return B::JFFltMat;
end
export Blmat!


end
