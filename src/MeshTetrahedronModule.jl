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
    return  T4block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt,:a)
end
export T4blocka

function T4blockb(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    return  T4block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt,:b)
end
export T4blockb

function T4blockca(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    return  T4block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt,:ca)
end
export T4blockca

function T4blockcb(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    return  T4block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt,:cb)
end
export T4blockcb

function T4block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt,orientation::Symbol)
    # Tetrahedral mesh of a rectangular block; orientation as given.
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
    return T4blockx(linspace(0.0,Length,nL+1),linspace(0.0,Width,nW+1),linspace(0.0,Height,nH+1),orientation);
end
export T4block

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

function  T4toT10(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetT4)
    # Convert a mesh of Tetrahedron T4 (four-node) to Tetrahedron T10.
    #
    # function [fens,fes] = T4_to_T10(fens,fes)
    #
    # Examples: 
    # [fens,fes] = T4_sphere(3.1,1);
    # [fens,fes] = T4_to_T10(fens,fes);
    # fens= onto_sphere(fens,3.1,connected_nodes(mesh_boundary(fes,[])));
    # figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
    #
    nedges=6;
    ec = [1  2; 2  3; 3  1; 4  1; 4  2; 4  3];
    # Additional node numbers are numbered from here
    newn=FENodeSetModule.count(fens)+1;
    # make a search structure for edges
    edges=MeshUtilModule.makecontainer();
    for i= 1:size(fes.conn,1)
        conn = fes.conn[i,:];
        for J = 1:nedges
            ev=conn[ec[J,:]];
            newn = MeshUtilModule.addhyperface!(edges, ev, newn);
        end
    end
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(JFFlt,newn-1,3);
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
            xyz[C[J].n,:]=mean(xyz[[i,C[J].o],:],1);
        end
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    # construct new geometry cells
    nconn=zeros(JFInt,size(fes.conn,1),10);
    nc=1;
    for i= 1:size(fes.conn,1)
        conn = fes.conn[i,:];
        econn=zeros(JFInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        nconn[nc,:] =[conn econn];
        nc= nc+ 1;
    end
    fes = FESetModule.FESetT10(conn=nconn);
    return fens,fes;
end
export T4toT10

function T10block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt; orientation::Symbol=:a)
# Tetrahedral (T10) mesh of a rectangular block.
#
# function [fens,fes] = T10_block(Length,Width,Height,nL,nW,nH)
#
# Range =<0,Length> x <0,Width> x <0,Height>
# Divided into elements: nL, nW, nH in the first, second, and
# third direction (x,y,z).
#
# See also: T4_blocka

    fens,fes = T4block(Length,Width,Height,nL,nW,nH,orientation);
    fens,fes = T4toT10(fens,fes);
    return fens,fes
end
export T10block

end
