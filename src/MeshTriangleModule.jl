module MeshTriangleModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule
using JFinEALE.MeshUtilModule

function T3blockx (xs::JFFltMat,ys::JFFltMat,orientation::Symbol=:a)
    return T3blockx (vec(xs),vec(ys),orientation)
end
export T3blockx

# T3 Mesh of a rectangle.
function T3blockx (xs::JFFltVec,ys::JFFltVec,orientation::Symbol=:a)
    #%
    #% Range =<0,Length> x <0,Width>
    #% Divided into triangular (t3) elements: nL, nW in the first, second (x,y).
    #% options = Attributes recognized by the constructor of fe_set_T3.
    #%
    #% Examples:
    #%     [fens,fes] = T3_block(3.1,2.2,3,2,[]);
    #%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
    #%
    #% See the visual gallery by running test_block
    #%
    #% See also:  T3_blocku, T3_cblock, T3_crossblock, T3_ablock,
    #%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn

    if (orientation!=:a) && (orientation!=:b)
        error("Cannot handle orientation :" * string(orientation))
    end
    
    nL =length(xs)-1;
    nW =length(ys)-1;
    nnodes=(nL+1)*(nW+1);
    ncells=2*(nL)*(nW);
    xys=zeros(JFFlt,nnodes,2);
    conns=zeros(JFInt,ncells,3);
    f=1;
    for j=1:(nW+1)
        for i=1:(nL+1)
            xys[f,:]=[xs[i] ys[j]];
            f=f+1;
        end
    end
    fens=FENodeSetModule.FENodeSet(xyz=xys);

    # function node_numbers1a(i,j,nL,nW)
    #     f=(j-1)*(nL+1)+i;
    #     nn=[f (f+1) f+(nL+1)];
    # end
    # function node_numbers2a(i,j,nL,nW)
    #     f=(j-1)*(nL+1)+i;
    #     nn=[(f+1)  f+(nL+1)+1 f+(nL+1)];
    # end
    # function nn=node_numbers1b(i,j,nL,nW) %#ok<INUSD>
    #     f=(j-1)*(nL+1)+i;
    #     nn=[f (f+1) f+(nL+1)+1];
    #     return;
    # end
    # function nn=node_numbers2b(i,j,nL,nW) %#ok<INUSD>
    #     f=(j-1)*(nL+1)+i;
    #     nn=[f  f+(nL+1)+1 f+(nL+1)];
    #     return;
    # end

    gc=1;
    for i=1:nL
        for j=1:nW
            f=(j-1)*(nL+1)+i;
            if     (orientation==:a)
                nn=[f (f+1) f+(nL+1)];
            elseif (orientation==:b)
                nn=[f (f+1) f+(nL+1)+1];
            end
            #nn=node_numbers1(i,j,nL,nW);
            conns[gc,:]=nn;
            gc=gc+1;
            if     (orientation==:a)
                nn=[(f+1)  f+(nL+1)+1 f+(nL+1)];
            elseif (orientation==:b)
                nn=[f  f+(nL+1)+1 f+(nL+1)];
            end
            #nn=node_numbers2(i,j,nL,nW);
            conns[gc,:]=nn;
            gc=gc+1;
        end
    end
    fes = FESetModule.FESetT3(conn=conns);
    
    return fens,fes
end
export T3blockx

# T3 Mesh of a rectangle.
function T3block (Length::JFFlt,Width::JFFlt,nL::JFInt,nW::JFInt,orientation::Symbol=:a)
    #%
    #% function [fens,fes] = T3_block(Length,Width,nL,nW,options)
    #%
    #% Range =<0,Length> x <0,Width>
    #% Divided into triangular (t3) elements: nL, nW in the first, second (x,y).
    #% options = Attributes recognized by the constructor of fe_set_T3.
    #%
    #% Examples:
    #%     [fens,fes] = T3_block(3.1,2.2,3,2,[]);
    #%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
    #%
    #% See the visual gallery by running test_block
    #%
    #% See also:  T3_blocku, T3_cblock, T3_crossblock, T3_ablock,
    #%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn

    
    return T3blockx (linspace(0.0,Length,nL+1),linspace(0.0,Width,nW+1),orientation)
end
export T3block

function T3toT6(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetT3)
# Convert a mesh of triangle T3 (three-node) to triangle T6.
#
# options =attributes recognized by the constructor fe_set_T6
    nedges=3;
    ec = [1 2; 2 3; 3 1];
    conns = fes.conn;
    # Additional node numbers are numbered from here
    newn=FENodeSetModule.count(fens)+1;
    # make a search structure for edges
    edges=MeshUtilModule.makecontainer();
    for i= 1:size(conns,1)
        conn = conns[i,:];
        for J = 1:nedges
            ev=conn[ec[J,:]];
            newn = MeshUtilModule.addhyperface!(edges, ev, newn);
        end
    end
     xyz1 =fens.xyz;             # Pre-existing nodes
   # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(JFFlt,newn-1,size(fens.xyz,2));
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
            xyz[C[J].n,:]=mean(xyz[[i,C[J].o],:],1);
        end
    end
   # construct new geometry cells
    nconns =zeros(JFInt,size(conns,1),6);
    nc=1;
    for i= 1:size(conns,1)
        conn = conns[i,:];
        econn=zeros(JFInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        nconns[nc,:] =[conn econn];
        nc= nc+ 1;
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    fes = FESetModule.FESetT6(conn=nconns) ;
    return fens,fes;
end
export T3toT6

function T6block(Length::JFFlt,Width::JFFlt,nL::JFInt,nW::JFInt,orientation::Symbol=:a)
    # Mesh of a rectangle of T6 elements. 
    # 
    # function [fens,fes] =T6_block(Length,Width,nL,nW,options)
    #  
    # Examples: 
    #     [fens,fes] = T6_block(3.1,2.2,3,2,[]);
    #     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
    #
    # If a different configuration of triangles is desired, run one of the
    # mesh generation functions below and then convert  the mesh to T6
    # triangles using T3_to_T6.
    #
    # See the visual gallery by running test_block
    # 
    # See also:  T3_blocku, T3_cblock,  T3_cblock2d, T3_crossblock, T3_ablock,
    #            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn, T3_to_T6

    fens,fes = T3block(Length,Width,nL,nW,orientation);
    fens,fes = T3toT6(fens,fes);
end
export T6block

end

