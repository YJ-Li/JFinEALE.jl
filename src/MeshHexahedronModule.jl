module MeshHexahedronModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule
using JFinEALE.MeshUtilModule
using JFinEALE.MeshModificationModule
using JFinEALE.MeshSelectionModule


function H8block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    # % Mesh of a 3-D block of H8 finite elements
    # %
    # % function [fens,fes] = H8_block(Length,Width,Height,nL,nW,nH)
    # %
    # % Arguments:
    # % Length,Width,Height= dimensions of the mesh in Cartesian coordinate axes,
    # % smallest coordinate in all three directions is  0 (origin)
    # % nL,nW,nH=number of elements in the three directions
    # %
    # % Range in xyz =<0,Length> x <0,Width> x <0,Height>
    # % Divided into elements: nL, nW, nH in the first, second, and
    # % third direction (x,y,z). Finite elements of type H8.
    # %
    # % Output:
    # % fens= finite element node set
    # % fes = finite element set
    # %
    # %
    # % Examples: 
    # %     [fens,fes] = H8_block(2,3,4, 2,6,4);
    # %     drawmesh({fens,fes},'fes','facecolor','red'); hold on
    # %
    # % See also: H8_block_u, H8_blockx, H8_composite_plate


    fens,fes = H8blockx(squeeze(linspace(0,Length,nL+1)',1),squeeze(linspace(0,Width,nW+1)',1),squeeze(linspace(0,Height,nH+1)',1));
end
export H8block

function H8blockx(xs::JFFltVec,ys::JFFltVec,zs::JFFltVec)
    # % Graded mesh of a 3-D block of H8 finite elements.
    #   %
    # % Mesh of a 3-D block, H8 cells. The nodes are located at the Cartesian
    # % product of the three intervals on the input.  This allows for
    # % construction of graded meshes.
    # %
    # % function [fens,fes] = H8_blockx(xs,ys,zs)
    # %
    # % xs,ys,zs =Locations of the individual planes of nodes. 
    # % 
    # % Output:
    # % fens= finite element node set
    # % fes = finite element set
    # %
    # %
    # % Examples:
    # % 
    # %  [fens,fes] = H8_blockx(1/125*(0:1:5).^3,4+(0:2:8),[3, 5, 9]); drawmesh ({fens,fes})
    # %  
    # %  [fens,fes] = H8_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2)
    # %  drawmesh({fens,fes},'fes','facecolor','red'); hold on
    # %
    # % See also: H8_block, H8_composite_plate

    nL =length(xs)-1;
    nW =length(ys)-1;
    nH =length(zs)-1;

    nnodes=(nL+1)*(nW+1)*(nH+1);
    xyz =zeros(JFFlt,nnodes,3);
    ncells=(nL)*(nW)*(nH);
    conns =zeros(JFInt,ncells,8);

    
    # first go along Length, then Width, and finally Height
    function node_numbers(i,j,k,nL,nW,nH)
        f=(k-1)*((nL+1)*(nW+1))+(j-1)*(nL+1)+i;
        nn=[f (f+1)  f+(nL+1)+1 f+(nL+1)];
        nn=[nn nn+((nL+1)*(nW+1))];
    end

    f=1;
    for k=1:(nH+1)
        for j=1:(nW+1)
            for i=1:(nL+1)
                xyz[f,:]=[xs[i] ys[j] zs[k]];
                f=f+1;
            end
        end
    end

    gc=1;
    for i=1:nL
        for j=1:nW
            for k=1:nH
                nn=node_numbers(i,j,k,nL,nW,nH);
                conns[gc,:]=nn;
                gc=gc+1;
            end
        end
    end
    # create the nodes
    fens = FENodeSetModule.FENodeSet(xyz= xyz);
    # Create the finite elements
    fes = FESetModule.FESetH8(conn=conns);
    
    return fens,fes;
end
export H8blockx


# % Create a solid mesh of 1/8 of the sphere of "radius". 
# %
# % function [fens,gcells]=H8_sphere(radius,nrefine)
# %
# % Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
# % four hexahedral elements if "nrefine==0", or more if "nrefine>0".
# % "nrefine" is the number of bisections applied  to refine the mesh.
# % 
function H8sphere(radius::JFFlt,nrefine::JFInt)
    a=sqrt(2.0)/2.0;
    b=1.0/sqrt(3.0);
    c=0.6*a;
    d=0.6*b;
    z= 0.0;
    h= 0.5;
    o= 1.0;
    xyz= [z  z  z;
        h  z  z;
        c  c  z;
        z  h  z;
        z  z  h;
        c  z  c;
        d  d  d;
        z  c  c;
        z  z  o;
        a  z  a;
        o  z  z;
        a  a  z;
        z  o  z;
        z  a  a;
        b  b  b]*radius;
    conns=[1  2  3  4  5  6  7  8;
        2  11  12  3  6  10  15  7;
        4  3  12  13  8  7  15  14;
        5  6  7  8  9  10  15  14];
  
    fens = FENodeSetModule.FENodeSet(xyz= xyz);
    fes = FESetModule.FESetH8(conn=conns);
    for i=1:nrefine
        fens,fes=H8refine(fens,fes);
        bg=MeshModificationModule.meshboundary(fes);
        l=MeshSelectionModule.feselect(fens,bg,facing=true, direction=[1,1,1]);
        cn = MeshSelectionModule.connectednodes(FESetModule.subset(bg,l))   ;
        for j=1:length(cn)
            fens.xyz[cn[j],:]=fens.xyz[cn[j],:]*radius/norm(fens.xyz[cn[j],:]);
        end
    end
     return fens, fes
end
export H8sphere

function H8refine(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetH8)
# % Refine a mesh of H8 hexahedrals by octasection.
# %
# % function [fens,fes] = H8_refine(fens,fes)
# %
# % Arguments and
# % Output:
# % fens= finite element node set
# % fes = finite element set
# %
# % Examples: 
# %
# %     xyz = [3, 1, 6; -5, 2, 1];
# %     [fens,fes] = H8_hexahedron(xyz,1,2,3);
# %     drawmesh({fens,fes},'fes','facecolor','red');
# %     [fens,fes] = H8_refine(fens,fes);
# %     figure;
# %     drawmesh({fens,fes},'fes','facecolor','m');

    fens,fes = H8toH27(fens,fes);
    conn=fes.conn;
    nconn=zeros(JFInt,8*size(conn,1),8);
    nc=1;
    for i= 1:size(conn,1)
        conn27=conn[i,:];
        nconn[nc,:] =conn27[[1,9,21,12,17,22,27,25]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[9,2,10,21,22,18,23,27]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[21,10,3,11,27,23,19,24]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[12,21,11,4,25,27,24,20]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[17,22,27,25,5,13,26,16]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[22,18,23,27,13,6,14,26]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[27,23,19,24,26,14,7,15]];        nc= nc+ 1;
        nconn[nc,:] =conn27[[25,27,24,20,16,26,15,8]];        nc= nc+ 1;
    end
   
    fes=FESetModule.FESetH8(conn=nconn);
    return fens, fes
end
export H8refine

function   H8toH27(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetH8)
# % Convert a mesh of hexahedra H8 to hexahedra H27.
# %
# % function [fens,fes] = H8_to_H27(fens,fes)
# %
# % Arguments and
# % Output:
# % fens= finite element node set
# % fes = finite element set
# %
    nedges=12;
    nfaces=6;
    ec = [1   2; 2   3; 3   4; 4   1; 5   6; 6   7; 7   8; 8   5; 1   5; 2   6; 3   7; 4   8;];
    fc = [1     4     3     2;
        1     2     6     5;
        2     3     7     6;
        3     4     8     7;
        4     1     5     8;
        6     7     8     5];
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
    # make a search structure for faces
    faces=MeshUtilModule.makecontainer();
    for i= 1:size(conns,1)
        conn = conns[i,:];
        for J = 1:nfaces
            fv=conn[fc[J,:]];
            newn = MeshUtilModule.addhyperface!(faces, fv, newn);
        end
    end
    # make a search structure for volumes
    volumes=MeshUtilModule.makecontainer();
    for i= 1:size(conns,1)
        conn = conns[i,:];
        newn = MeshUtilModule.addhyperface!(volumes, conn, newn);
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
   # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(faces)
        C=faces[i];
        for J = 1:length(C)
            xyz[C[J].n,:] = mean(xyz[[i,C[J].o],:],1);
        end
    end
    # calculate the locations of the new nodes
    # and construct the new nodes
    for i in keys(volumes)
       C=volumes[i];
        for J = 1:length(C)
            xyz[C[J].n,:] = mean(xyz[[i,C[J].o],:],1);
        end
    end
     # construct new geometry cells
    nconns =zeros(JFInt,size(conns,1),27);
    nc=1;
    for i= 1:size(conns,1)
        conn = conns[i,:];
        econn=zeros(JFInt,1,nedges);
        for J = 1:nedges
            ev=conn[ec[J,:]];
            h,n=MeshUtilModule.findhyperface!(edges, ev);
            econn[J]=n;
        end
        fconn=zeros(JFInt,1,nfaces);
        for J = 1:nfaces
            fv=conn[fc[J,:]];
            h,n=MeshUtilModule.findhyperface!(faces, fv);
            fconn[J]=n;
        end
        h,n=MeshUtilModule.findhyperface!(volumes, conn);
        vconn=n;
        nconns[nc,:] =[conn econn fconn vconn];
        nc= nc+ 1;
    end
    fens =FENodeSetModule.FENodeSet(xyz);
    fes = FESetModule.FESetH27(conn=nconns) ;
    return fens,fes;
end
export H8toH27

function H8hexahedron(xyz::JFFltMat,nL::JFInt,nW::JFInt,nH::JFInt;block_mesh_handle=nothing)
# % Mesh of a general hexahedron given by the location of the vertices.
# %
# % function [fens,fes] = H8_hexahedron(xyz,nL,nW,nH,block_mesh_handle)
# % 
# % xyz = One vertex location per row; Either two rows (for a rectangular
# %      block given by the its corners), or eight rows (general hexahedron).
# % nL, nW, nH = Divided into elements: nL, nW, nH in the first, second, and
# %      third direction. 
# % Optional argument:
# % block_mesh_handle = function handle of the block-generating mesh function
# %      (having the signature of the function H8_block()).
# %
# % Output:
# % fens= finite element node set
# % fes = finite element set
# %
# %
# % Examples: 
# %
# %     xyz = [3, 1, 6; -5, 2, 1];
# %     [fens,fes] = H8_hexahedron(xyz,12,3,4);
# %     drawmesh({fens,fes},'fes','facecolor','red'); hold on
# %
# %     A=[0,0,0]; B=[0,0,2]; C=[0,3,2]; D=[0,3,0];
# %     E=[5,0,0]; F=[5,0,2]; G=[5,3,2]; H=[5,3,0];
# %     P=[3.75,0,0];
# %     [fens,fes] = H8_hexahedron([A;P;(D+H)/2;D;B;(B+F)/2;(C+G)/2;C],2,3,4,[]);
# %     drawmesh({fens,fes},'fes','facecolor','red'); hold on
# %
# %     A=[0,0,0]; B=[0,0,2]; C=[0,3,2]; D=[0,3,0];
# %     E=[5,0,0]; F=[5,0,2]; G=[5,3,2]; H=[5,3,0];
# %     P=[3.75,0,0];
# %     [fens,fes] = H8_hexahedron([A;P;(D+H)/2;D;B;(B+F)/2;(C+G)/2;C],1,2,3,@H20_block);
# %     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on

    npts=size(xyz,1);
    if npts==2
        lo=minimum(xyz,1);
        hi=maximum(xyz,1);
        xyz=[lo[1]  lo[2]  lo[3];
            hi[1]  lo[2]  lo[3];
            hi[1]  hi[2]  lo[3];
            lo[1]  hi[2]  lo[3];
            lo[1]  lo[2]  hi[3];
            hi[1]  lo[2]  hi[3];
            hi[1]  hi[2]  hi[3];
            lo[1]  hi[2]  hi[3]];
    elseif npts!=8
        error("Need 2 or 8 points");
    end

    if block_mesh_handle==nothing
        block_mesh_handle =H8block; # default block type
    end

    fens,fes= block_mesh_handle(2.0,2.0,2.0,nL,nW,nH);

    dummy=FESetModule.FESetH8(conn=reshape([1:8],1,8))
    pxyz=fens.xyz;
    for i=1:FENodeSetModule.count(fens)
        N = FESetModule.bfun(dummy,pxyz[i,:]-1.0);# shift coordinates by -1
        pxyz[i,:] =N'*xyz;
    end
    fens.xyz=pxyz;
    
    return fens,fes;
end
export H8hexahedron

function H27block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
    # % Mesh of a 3-D block of H27 finite elements
    # %
    # % Arguments:
    # % Length,Width,Height= dimensions of the mesh in Cartesian coordinate axes,
    # % smallest coordinate in all three directions is  0 (origin)
    # % nL,nW,nH=number of elements in the three directions
    # %
    # % Range in xyz =<0,Length> x <0,Width> x <0,Height>
    # % Divided into elements: nL, nW, nH in the first, second, and
    # % third direction (x,y,z). Finite elements of type H27.
    # %
    # % Output:
    # % fens= finite element node set
    # % fes = finite element set
    # %    
    # %
    # % Examples: 
    # %
    # %     [fens,fes] = H27_block(2,3,4, 1,2,3);
    # %     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
    # %
    # % See also: H8_block, H8_to_H27
    # %
    fens,fes = H8block(Length,Width,Height,nL,nW,nH);
    fens,fes = H8toH27(fens,fes);
    return fens,fes
end
export H27block

function doextrude(fens,fes::FESetQ4,nLayers,extrusionh)
    nn1=count(fens);
    nnt=nn1*nLayers;
    ngc=count(fes)*nLayers;
    hconn=zeros(JFInt,ngc,8);
    xyz =zeros(JFFlt,nn1*(nLayers+1),3);
    for j=1:nn1
        xyz[j,:] =extrusionh(fens.xyz[j,:],0);
    end
    for k=1:nLayers
        for j=1:nn1
            f=j+k*nn1;
            xyz[f,:] =extrusionh(fens.xyz[j,:],k);
        end
    end
    
    gc=1;
    for k=1:nLayers
        for i=1:count(fes)
            hconn[gc,:]=[fes.conn[i,:]+(k-1)*nn1 fes.conn[i,:]+k*nn1];
            gc=gc+1;
        end
    end
    efes = FESetModule.FESetH8(conn=hconn);
    efens = FENodeSetModule.FENodeSet(xyz= xyz);
    return efens,efes
end

function H8extrudeQ4(fens::FENodeSet, fes::FESetQ4,nLayers::JFInt,extrusionh::Function)
    # % Extrude a mesh of quadrilaterals into a mesh of hexahedra (H8).
    # %
    # % function [fens,fes] = H8_extrude_Q4(fens,fes,nLayers,extrusionh)
    #     %
    #     % Inputs:
    #     % fens,fes = nodes and cells of the quadrilateral mesh,
    #     % nLayers= number of layers,
    #     % extrusionh= handle to function that computes the extrusion vector
    #         %   nxyz=extrusionh(xyz,k)
    #         % 
    #         % Output:
    #         % fens= node set which incorporates the nodes on input plus the nodes of
    #         % the extruded elements. Note that no merging (fusing) of nodes is performed.
    #         % h8fes = H8 extruded geometric cells
    #         %
    #         % Examples: 
    #         %     L= 0.3; % in-plane dimension
    #         %     W = 0.3; % in-plane dimension
    #         %     a= 0.15; % hole radius
    #         %     H = 0.03; % thickness of the plate
    #         %     nL=5;nH=1;nW=5;na=7;
    #         %     [fens,fes]=Q4_elliphole(a,a,L,W,nL,nW,na,[]);
    #         %     [fens,fes] = H8_extrude_Q4(fens,fes,nH,@(x,i)([x,0]+[0,0,H*i]));
    #         %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    #         %
    id = [1:count(fens)];
    cn=connectednodes(fes);
    id[cn[:]]=vec([1:length(cn)]);
    q4fes= deepcopy(fes);
    updateconn!(q4fes,id);
    q4fens = FENodeSetModule.FENodeSet(xyz= fens.xyz[cn[:],:]);
    h8fens,h8fes= doextrude(q4fens,q4fes,nLayers,extrusionh);
    return h8fens,h8fes
end
export H8extrudeQ4

function H8spheren(radius::JFFlt,nperradius::JFInt)
# % Create a solid mesh of 1/8 of sphere. 
# %
# % Create a solid mesh of 1/8 of the sphere of "radius", with nperradius 
# % elements per radius. 
# %
# % function [fens,fes]=H8_sphere_n(radius,nperradius)
# %
# % Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
# % 4*(nperradius/2)^2 hexahedral elements.
# % 
# % Output:
# % fens= finite element node set
# % fes = finite element set
# %
# % Examples: 
# %     [fens,fes]=H8_sphere_n(22.3,3);
# %     drawmesh({fens,fes},'fes','facecolor','red'); hold on
# %
# % See also: H8_sphere


    if (mod(nperradius,2)!=0)
        nperradius=nperradius+1;
    end
    nL=convert(JFInt,nperradius/2); nW=nL; nH=nL;
  
    a=sqrt(2.0)/2.0;
    b=1.0/sqrt(3.0);
    c=0.6*a;
    d=0.6*b;
    z= 0.0;
    h= 0.5;
    o= 1.0;
    xyz= [z  z  z;
        h  z  z;
        c  c  z;
        z  h  z;
        z  z  h;
        c  z  c;
        d  d  d;
        z  c  c;
        z  z  o;
        a  z  a;
        o  z  z;
        a  a  z;
        z  o  z;
        z  a  a;
        b  b  b]*radius;
    conns=[1  2  3  4  5  6  7  8;
        2  11  12  3  6  10  15  7;
        4  3  12  13  8  7  15  14;
           5  6  7  8  9  10  15  14];
    tolerance=radius*1.0e-6;
    
  # fens = FENodeSetModule.FENodeSet(xyz= xyz);
  #   fes = FESetModule.FESetH8(conn=conns);
    fens,fes = H8hexahedron(xyz[conns[1,:][:],:],nL,nW,nH);
    fens1,fes1 = H8hexahedron(xyz[conns[2,:][:],:],nL,nW,nH);
    fens,fes1,fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance);
    fes=cat(fes1,fes2);
    fens1,fes1 = H8hexahedron(xyz[conns[3,:][:],:],nL,nW,nH);
    fens,fes1,fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance);
    fes=cat(fes1,fes2);
    fens1,fes1 = H8hexahedron(xyz[conns[4,:][:],:],nL,nW,nH);
   fens,fes1,fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance);
    fes=cat(fes1,fes2);
    
    xyz=deepcopy(fens.xyz);
    layer=oftype(1.0,Inf) + zeros(JFFlt,size(xyz, 1),1);
    conn=deepcopy(fes.conn);
    bg=meshboundary(fes);
    l=feselect(fens,bg;facing=true, direction=[1. 1. 1.]);
    cn = connectednodes(subset(bg,l))   ;
    layer[cn]=1;
    for j=1:nperradius-1
        for k=1:size(conn,1)
            ll=layer[conn[k,:]];
            ml=minimum(ll);
            if (ml==j)
                ix=isinf(ll);
                ll[ix]=j+1;
                layer[conn[k,:]]=ll;
            end
        end 
    end
    nxyz=deepcopy(xyz);
    for j=1:size(xyz,1)
        if (!isinf(layer[j]))
            nxyz[j,:]=nxyz[j,:]*(nperradius-layer[j]+1)/nperradius*radius/norm(nxyz[j,:]);
        end
    end
    s= linspace(0., 1., nperradius);
    for j=1:size(xyz,1)
        if (!isinf(layer[j]))
            nxyz[j,:]=s[layer[j]]*xyz[j,:]+(1-s[layer[j]])*nxyz[j,:];
        end
    end
    fens.xyz=deepcopy(nxyz);
    return fens,fes
end
export H8spheren

function H20block(Length::JFFlt,Width::JFFlt,Height::JFFlt,nL::JFInt,nW::JFInt,nH::JFInt)
# Mesh of a 3-D block of H20 finite elements
#
# Arguments:
# Length,Width,Height= dimensions of the mesh in Cartesian coordinate axes,
# smallest coordinate in all three directions is  0 (origin)
# nL,nW,nH=number of elements in the three directions
#
# Range in xyz =<0,Length> x <0,Width> x <0,Height>
# Divided into elements: nL, nW, nH in the first, second, and
# third direction (x,y,z). Finite elements of type H20.
#
# Output:
# fens= finite element node set
# fes = finite element set
#
#
# Examples: 
#     [fens,fes] = H20_block(2,3,4, 1,2,3);
#     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
#
# See also: H8_block, H8_to_H20
#
    fens,fes = H8block(Length,Width,Height,nL,nW,nH);
    fens,fes = H8toH20(fens,fes);
end
export H20block

function   H8toH20(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESetH8)
# % Convert a mesh of hexahedra H8 to hexahedra H20.
# %
# % Arguments and
# % Output:
# % fens= finite element node set
# % fes = finite element set
# %
    nedges=12;
    ec = [1   2; 2   3; 3   4; 4   1; 5   6; 6   7; 7   8; 8   5; 1   5; 2   6; 3   7; 4   8;];
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
   # construct new geometry cells
    nconns =zeros(JFInt,size(conns,1),20);
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
    fes = FESetModule.FESetH20(conn=nconns) ;
    return fens,fes;
end
export H8toH20

end

