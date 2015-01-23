module MeshTetrahedronModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule
using JFinEALE.MeshUtilModule
using JFinEALE.MeshModificationModule
using JFinEALE.MeshSelectionModule

# function T4cylinderdel(Length::JFFlt,Radius::JFFlt,nL::JFInt,nR::JFInt)
# # Tetrahedral (T4) Delaunay Mesh of a cylinder.
# #
# # function [fens,fes] = T4_cylinderdel(Length,Radius,nL,nR)
# #
# # Tetrahedral (T4) Delaunay Mesh of a cylinder of Length and Radius.
# # Divided into elements: nL (lengthwise), nR (along the radius).
# #
# # Examples: 
# # [fens,fes] = T4_cylinderdel(3.1,1.7,5,2);
# # figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on

#     dR=Radius/(nR+1/2);
#     xnnodes= 3*(2*(0:1:nR)+1);
#     nnodes=(nL+1)*sum(xnnodes);
#     x=zeros(nnodes,1);
#     y=zeros(nnodes,1);
#     z=zeros(nnodes,1);
#     xs=zeros(nnodes,3);
#     k = 1;
#     for j = 0:1:nL
#         xjiggle = Length/nL/10*j*(j-nL)/nL/nL*rand;
#         for n= 0:1:nR
#             r =dR*(1/2+n);
#             rjiggle=Radius/nR/10*rand*(nR-n)/nR;
#             r= r+rjiggle;
#             dA =2*pi/xnnodes(n+1);
#             dAjiggle =dA/10*rand*(nR-n)/nR;
#             for m=1:xnnodes(n+1)
#                 x(k)=j*Length/nL+xjiggle;
#                 y(k)=r*cos(dAjiggle+dA*(m-1));
#                 z(k)=r*sin(dAjiggle+dA*(m-1));
#                 xs(k,:)=[x(k),y(k),z(k)];
#                 k=k+1;
#             end
#         end
#     end
#     if (verLessThan('matlab', '7.6'))
#         T = delaunay3(x,y,z);# 
#     else
#         T = DelaunayTri(x,y,z);# replaced DELAUNAY3 that will be removed in a future release
#     end

#     ncells=size(T,1);
#     conns=zeros(ncells,4);
#     for i=1:ncells
#         xyz=zeros (4, 3);
#         for k=1:4
#             xyz(k,:) = xs(T(i,k),:);
#         end
#         if tetvol(xyz) > 0
#             conns(i,:) =T(i,:);
#         else
#             conns(i,:) =T(i,[1, 3, 2, 4]);
#         end
#     end
#     fens=fenode_set(struct('xyz',xs));
#     fes =fe_set_T4(struct ('conn',conns));
# end

# #compute volume of a tetrahedron
# # Given the 4x3 vertex coordinate matrix V of a tetrahedron, TETVOL(V)
# # returns the volume of the tetrahedron.
# function vol = tetvol(v)
#     vol = det([v(2,:)-v(1,:);v(3,:)-v(1,:);v(4,:)-v(1,:)])/6;
#     #     if abs (vol) < 0.1
#     #         warning (' sliver?')
#     #     end
#     return;
# end

function T4blocka(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    # Tetrahedral mesh of a rectangular block; orientation 'a'.
    #
    # function [fens,fes] = T4_blocka(Length,Width,Height,nL,nW,nH)
    #
    # The mesh is produced by splitting each logical  rectangular cell into six
    # tetrahedra.
    # Range =<0,Length> x <0,Width> x <0,Height>
    # Divided into elements: nL, nW, nH in the first, second, and
    # third direction (x,y,z).
    #
    # Output:
    # fens= finite element node set
    # fes = finite element set
    #
    # Examples:
    #         [fens,fes] = T4_blocka(2,3,4, 2,6,4);
    #         figure; drawmesh({fens,fes},'fes','facecolor','red'); hold on
    #
    # See also: T4_blockx,  T4_blockb,  T4_blockca,  T4_blockcb
    return T4blockx(linspace(0.0,Length,nL+1),linspace(0.0,Width,nW+1),linspace(0.0,Height,nH+1),:a);
end
export T4blocka

function T4blockx(xs::JFFltMat,ys::JFFltMat,zs::JFFltMat,orientation::Symbol)
    return T4blockx(vec(xs),vec(ys),vec(zs),orientation)
end

function T4blockx(xs::JFFltVec,ys::JFFltVec,zs::JFFltVec,orientation::Symbol)
    # Graded tetrahedral (T4) mesh of a rectangular block.
    #
    # function [fens,fes] = T4_blockx(xs,ys,zs,orientation)
    #
    # xs,ys,zs =Locations of the individual planes of nodes.
    # orientation = either one of 'a', 'b', 'ca', 'cb' (when omitted, default is 'a')
    #
    # The mesh is produced by splitting each logical  rectangular cell into six
    # tetrahedra.
    #
    # Output:
    # fens= finite element node set
    # fes = finite element set
    #
    # Examples:
    #     [fens,fes] = T4_blockx(1/125*(0:1:5).^3,4+(0:2:8),[3, 5, 9]);
    #     drawmesh ({fens,fes})
    #
    #     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'a');
    #     figure; drawmesh({fens,fes},'fes','facecolor','red'); hold on
    #
    #     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'b');
    #     figure; drawmesh({fens,fes},'fes','facecolor','g'); hold on
    #
    #     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'ca');
    #     figure; drawmesh({fens,fes},'fes','facecolor','b'); hold on
    #
    #     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'cb');
    #     figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
    #
    # See also: T4_blocka,  T4_blockb,  T4_blockca,  T4_blockcb
    nL =length(xs)-1;
    nW =length(ys)-1;
    nH =length(zs)-1;
    nnodes=(nL+1)*(nW+1)*(nH+1);
    ncells=6*(nL)*(nW)*(nH);
    xyzs=zeros(JFFlt,nnodes,3);
    conns=zeros(JFInt,ncells,4);
    if (orientation==:a)
        t4ia = [1  8  5  6; 3  4  2  7; 7  2  6  8; 4  7  8  2; 2  1  6  8; 4  8  1  2];
        t4ib = [1  8  5  6; 3  4  2  7; 7  2  6  8; 4  7  8  2; 2  1  6  8; 4  8  1  2];
    elseif (orientation==:b)
        t4ia = [2 7 5 6; 1 8 5 7; 1 3 4 8; 2 1 5 7; 1 2 3 7; 3 7 8 1];
        t4ib = [2 7 5 6; 1 8 5 7; 1 3 4 8; 2 1 5 7; 1 2 3 7; 3 7 8 1];
    elseif (orientation==:ca)
        t4ia = [8  4  7  5; 6  7  2  5; 3  4  2  7; 1  2  4  5; 7  4  2  5];
        t4ib = [7  3  6  8; 5  8  6  1; 2  3  1  6; 4  1  3  8; 6  3  1  8];
    elseif (orientation==:cb)
        t4ia = [7  3  6  8; 5  8  6  1; 2  3  1  6; 4  1  3  8; 6  3  1  8];
        t4ib = [8  4  7  5; 6  7  2  5; 3  4  2  7; 1  2  4  5; 7  4  2  5];
    end
    f=1;
    for k=1:(nH+1)
        for j=1:(nW+1)
            for i=1:(nL+1)
                xyzs[f,:]=[xs[i] ys[j] zs[k]];
                f=f+1;
            end
        end
    end

    fens=FENodeSetModule.FENodeSet(xyz= xyzs);
    
    function node_numbers(i::JFInt,j::JFInt,k::JFInt,nL::JFInt,nW::JFInt,nH::JFInt)
        f=(k-1)*((nL+1)*(nW+1))+(j-1)*(nL+1)+i;
        nn=[f (f+1)  f+(nL+1)+1 f+(nL+1)];
        return [nn nn+((nL+1)*(nW+1))];
    end

    gc=1;
    for i=1:nL
        for j=1:nW
            for k=1:nH
                nn=node_numbers(i,j,k,nL,nW,nH);
                if (mod (sum( [i,j,k] ),2)==0)
                    t4i =t4ib;
                else
                    t4i =t4ia;
                end
                for r=1:size(t4i,1)
                    for c1=1:size(t4i,2)
                        conns[gc,c1]=nn[t4i[r,c1]];
                    end                    
                    gc=gc+1;
                end
            end
        end
    end
    fes = FESetModule.FESetT4(conn=conns[1:gc-1,:]);
    
    return fens,fes
end
export T4blockx

end
