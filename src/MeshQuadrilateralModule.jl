module MeshQuadrilateralModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule
using JFinEALE.MeshModificationModule
using JFinEALE.MeshUtilModule


function Q4annulus(rin::JFFlt,rex::JFFlt,nr::JFInt,nc::JFInt,Angl::JFFlt)
    # % Mesh of an annulus segment.
    # %
    # % function [fens,fes] = Q4_annulus(rin,rex,nr,nc,thickness)
    # %
    # % Mesh of an annulus segment, centered at the origin, with internal radius
    # % rin, and  external radius rex, and  development angle Angl. Divided into
    # % elements: nr, nc in the radial and circumferential direction
    # % respectively.
    # %
    # % Note that if you wish to have an annular region with 360° development
    # % angle  (closed annulus), the nodes along the slit  need to be merged.
    # %
    # % Examples: 
    # %     [fens,fes] = Q4_annulus(1.5,2.5,2,17,1.0);
    # %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    # %
    # % See also: Q4_block
    # %
    trin=min(rin,rex);
    trex=max(rin,rex);
    fens,fes =Q4block(trex-trin,Angl,nr,nc);
    xy=fens.xyz;
    for i=1:FENodeSetModule.count(fens)
        r=trin+xy[i,1]; a=xy[i,2];
        xy[i,:]=[r*cos(a) r*sin(a)];
    end
    fens.xyz=xy;
    return fens,fes
end
export Q4annulus

function Q4quadrilateral(xyz::JFFltMat, nL::JFInt, nW::JFInt)
    
    # % Mesh of a general quadrilateral given by the location of the vertices.
    # %
    # % function [fens,fes] = Q4_quadrilateral(xyz,nL,nW,options)
    # %
    # % xyz = One vertex location per row; Either two rows (for a rectangular
    # % block given by the two corners), or four rows (General quadrilateral).
    # % Divided into elements: nL, nW in the first and second direction.
    # % options = Attributes recognized by the constructor of fe_set_Q4.
    # %
    # % Examples: 
    # % [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
    # % drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    # %
    npts=size(xyz,1);
    if npts==2
        lo=minimum(xyz,1);
        hi=maximum(xyz,1);
        xyz=[[lo[1] lo[2]];
             [hi[1] lo[2]];
             [hi[1] hi[2]];
             [lo[1] hi[2]]];
    elseif npts!=4
        error("Need 2 or 4 points");
    end

    fens,fes = Q4block(2.,2.,nL,nW);

    xyz1=fens.xyz;
    if (size(xyz1,2)<size(xyz,2))
        nxyz1=zeros(JFFlt,size(xyz1,1),size(xyz,2));
        nxyz1[:,1:size(xyz1,2)]=xyz1;
        xyz1=nxyz1;
    end
    
    dummy=FESetModule.FESetQ4(conn=reshape([1:4],1,4))
    pxyz=fens.xyz;
    for i=1:FENodeSetModule.count(fens)
        N = FESetModule.bfun(dummy,pxyz[i,:]-1.0);# shift coordinates by -1
        pxyz[i,:] =N'*xyz;
    end
    fens.xyz=xyz1;
    return fens,fes 
end
export Q4quadrilateral

function Q4elliphole(xradius::JFFlt,yradius::JFFlt,L::JFFlt,H::JFFlt,nL::JFInt,nH::JFInt,nR::JFInt)
    # % Mesh of one quarter of a rectangular plate with an elliptical hole
    # %
    # % function [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options)
    # %
    # % xradius,yradius = radius of the ellipse,
    # % L,H= and dimensions of the plate,
    # % nL,nH= numbers of edges along the side of the plate,
    # % nR= number of edges along the circumference,
    # % options= options accepted by fe_set_Q4
    # %
    # % Examples: 
    # %     [fens,fes]=Q4_elliphole(1.2,2.4,4.8,3.5,4,2,2,[]);
    # %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    # %
    # %     [fens,fes]=Q4_elliphole(2.4,1.2,4.8,3.5,4,2,6,[]);
    # %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on

    dA =pi/2/(nL +nH);
    tolerance =(xradius+yradius)/(nL*nH)/100;
    fens= nothing; fes= nothing;
    for i= 1:nH
        xy = [xradius*cos((i-1)*dA) yradius*sin((i-1)*dA);
              L (i-1)/nH*H;
              L (i)/nH*H;
              xradius*cos((i)*dA) yradius*sin((i)*dA)];
        fens1,fes1 = Q4quadrilateral(xy,nR,1);
        if (fens== nothing)
            fens=fens1; fes =fes1;
        else
            fens,fes1,fes2 = MeshModificationModule.mergemeshes(fens1, fes1, fens, fes, tolerance);
            fes =FESetModule.cat(fes1,fes2);
        end
    end
    for i= 1:nL
        xy = [xradius*cos((nH+i-1)*dA)   yradius*sin((nH+i-1)*dA);
              (nL-i+1)/nL*L   H;
              (nL-i)/nL*L  H;
              xradius*cos((nH+i)*dA)   yradius*sin((nH+i)*dA)];
        fens1,fes1 = Q4quadrilateral(xy,nR,1);
        fens,fes1,fes2 = MeshModificationModule.mergemeshes(fens1, fes1, fens, fes, tolerance);
        fes =FESetModule.cat(fes1,fes2);
    end
    return fens,fes
end
export Q4elliphole

function Q4block(Length::JFFlt,Width::JFFlt,nL::JFInt,nW::JFInt)
    # % Mesh of a rectangle, Q4 elements
    # %
    # % function [fens,fes] = Q4_block(Length,Width,nL,nW,options)
    # %
    # % Rectangle <0,Length> x <0,Width>
    # % Divided into elements: nL, nW in the first, second (x,y).
    # % options = structure with fields recognized by the constructor of the
    # %       Q4 finite element
    # % 
    # % Examples: 
    # % [fens,fes] = Q4_block(3.5,1.75,2,3,struct('other_dimension',1))
    # % drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    # %  
    # % See also: Q4_blockx, fe_set_Q4

    fens,fes = Q4blockx(squeeze(linspace(0,Length,nL+1)',1),squeeze(linspace(0,Width,nW+1)',1));
end
export Q4block

function Q4blockx(xs::JFFltVec,ys::JFFltVec)
    # % Graded mesh  of a rectangle, Q4 finite elements.
    # %
    # % function [fens,fes] = Q4_blockx(xs, ys, options)
    # %
    # % Mesh of a 2-D block, Q4 finite elements. The nodes are located at the
    # % Cartesian product of the two intervals on the input.  This allows for
    # % construction of graded meshes.
    # %
    # % xs,ys - Locations of the individual planes of nodes.
    # %
    # % options - structure with fields recognized by the constructor of the
    # %   fe_set_Q4 object
    # % 
    # % Examples:  
    # %     [fens,fes] = Q4_blockx(1/125*(0:1:7).^3,4+(0:2:8), []);
    # %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    # %  
    # % See also: Q4_block, fe_set_Q4
    # %

    nL = length(xs) - 1;
    nW = length(ys) - 1;

    nnodes = (nL+1) * (nW+1);
    ncells = nL * nW;

    # preallocate node locations
    xyz = zeros(JFFlt, nnodes, 2);
    k = 1;
    for j = 1:(nW+1)
        for i = 1:(nL+1)
            xyz[k,:] = [xs[i] ys[j]];
            k = k + 1;
        end
    end
    # create the nodes
    fens = FENodeSetModule.FENodeSet(xyz= xyz);

    #preallocate connectivity matrix
    conn = zeros(JFInt, ncells, 4);
    
    function  nodenumbers(i,j,nL,nW)
        f = (j-1) * (nL+1) + i;
        nn = [f, (f+1), f+(nL+1)+1, f+(nL+1)];
        return nn
    end
    
    k = 1;
    for i = 1:nL
        for j = 1:nW
            conn[k,:] = nodenumbers(i,j,nL,nW);
            k = k + 1;
        end
    end
    # create the cells
    fes = FESetModule.FESetQ4(conn=conn);
    
    return fens,fes;
end
export Q4blockx

function Q8block(Length::JFFlt,Width::JFFlt,nL::JFInt,nW::JFInt)
    # Mesh of a rectangle of Q8 elements.
    # 
    # function [fens,fes] =Q8_block(Length,Width,nL,nW,options)
    # 
    # Examples: 
    # [fens,fes] = Q8_block(3.5,1.75,2,3,struct('other_dimension',1))
    # drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    #
    # See also: Q4_block, Q4_to_Q8

    fens,fes  = Q4block(Length,Width,nL,nW);
    fens,fes = Q4toQ8(fens,fes);
end
export Q8block

function Q4toQ8(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetQ4)
    # Convert a mesh of quadrilateral Q4 to quadrilateral Q8.
    #
    # function [fens,fes] = Q4_to_Q8(fens,fes,options)
    #
    # options =attributes recognized by the constructor fe_set_Q8
    #
    # Examples: 
    #     R=8.5;
    #     [fens,fes]=Q4_sphere(R,1,1.0);
    #     [fens,fes] = Q4_to_Q8(fens,fes,[]);
    #     fens= onto_sphere(fens,R,[]);
    #     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    #
    # See also: fe_set_Q8

    nedges=4;
    ec = [1  2; 2  3; 3  4; 4  1];
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
    xyz =zeros(JFFlt,newn-1,size(xyz1,2));
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
    nconns =zeros(JFInt,size(conns,1),8);
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
    fes = FESetModule.FESetQ8(conn=nconns);
    return fens,fes
end
export Q4toQ8


# Refine a mesh of quadrilaterals by bisection
#
# function [fens,fes] = Q4_refine(fens,fes)
#
# Examples: 
# [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
# [fens,fes] = Q4_refine(fens,fes);
# drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on

function Q4refine(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetQ4)
    nedges=4;
    ec = [1  2; 2  3; 3  4; 4  1];
    # make a search structure for edges
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
    newn=  newn+size(fes.conn,1) # add the interior nodes to the total
    xyz1 =fens.xyz;             # Pre-existing nodes
    # Allocate for vertex nodes plus edge nodes plus face nodes
    xyz =zeros(JFFlt,newn-1,size(xyz1,2));
    xyz[1:size(xyz1,1),:] = xyz1; # existing nodes are copied over
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(edges)
        C=edges[i];
        for J = 1:length(C)
            xyz[C[J].n,:]=mean(xyz[[i,C[J].o],:],1);
        end
    end
    # construct new geometry cells: for new elements out of one old one
    nconn =zeros(JFInt,4*size(fes.conn,1),4);
    nc=1;
    for i= 1:size(fes.conn,1)
        conn = fes.conn[i,:];
        econn=zeros(JFInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        inn=size(xyz,1)-size(fes.conn,1)+i
        xyz[inn,:]=mean(xyz[conn[:],:],1); # interior node
        #h,inn=MeshUtilModule.findhyperface!(faces, conn);
        nconn[nc,:] =[conn[1] econn[1] inn econn[4]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[2] econn[2] inn econn[1]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[3] econn[3] inn econn[2]];
        nc= nc+ 1;
        nconn[nc,:] =[conn[4] econn[4] inn econn[3]];
        nc= nc+ 1;
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    nfes = FESetModule.FESetQ4(conn=nconn);
    return fens,nfes            # I think I should not be overwriting the input!
end
export Q4refine

# Convert a mesh of quadrilateral Q4's to two T3 triangles  each.
#
# function [fens,fes] = Q4_to_T3(fens,fes,options)
#
# options =attributes recognized by the constructor fe_T3, and
# orientation = 'default' or 'alternate' chooses which diagonal is taken
#      for splitting
# Example: 
#     [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
#     [fens,fes] = Q4_to_T3(fens,fes,[]);
#     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
#
#     [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
#     [fens,fes] = Q4_to_T3(fens,fes,struct('orientation','alternate'));
#     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
# will generate triangles by using the alternate diagonal for splitting.
# 
# See also: Q4_to_T3_sd

function Q4toT3(fens,fes::FESetQ4, orientation:: Symbol=:default)
    connl1=[1, 2, 3];
    connl2=[1, 3, 4];
    if orientation==:alternate
        connl1=[1, 2, 4];
        connl2=[3, 4, 2];
    end
    nconns=zeros(JFInt,2*count(fes),3);
    nc=1;
    for i= 1:count(fes)
        conn = fes.conn[i,:];
        nconns[nc,:] =conn[connl1];
        nc= nc+ 1;
        nconns[nc,:] =conn[connl2];
        nc= nc+ 1;
    end
    nfes = FESetModule.FESetT3(conn=nconns);
    return fens,nfes            # I think I should not be overwriting the input!
end
export Q4toT3

end
