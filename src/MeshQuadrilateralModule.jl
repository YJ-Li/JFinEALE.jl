module MeshQuadrilateralModule

using JFFoundationModule
using FESetModule
using FENodeSetModule
using MeshModificationModule


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


end

