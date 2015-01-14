module MeshLineModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule
using JFinEALE.MeshModificationModule


function L2block(Length::JFFlt,nL::JFInt)
    # % Mesh of a 1-D block of H8 finite elements
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


    fens,fes = L2blockx(squeeze(linspace(0,Length,nL+1)',1));
end
export L2block

function L2blockx(xs::JFFltVec)
# % Graded mesh of a 1-D block, L2 finite elements. 
# %
# % Syntax
# %
# % function [fens,fes] = L2_blockx(xyz,options)
# %
# % Input
# %
# % xyz - (n x p) matrix containing the locations of the individual nodes.
# %
# % options - structure with fields recognized by the constructor of the
# %   L2 finite element
# %
# % Example
# % 
# %     [fens,fes] = L2_blockx((2:1:7)'.^3, struct('other_dimension', 0.1));
# %     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
# %
# %     a=pi/180*(0:1:6)'/6*135;
# %     [fens,fes] = L2_blockx(5.1*[sin(a),cos(a)], struct('other_dimension', 0.1));
# %     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
# % See also: L2_block, fe_set_L2

    xyz =reshape(sort(xs),length(xs),1);
    ncells=length(xs)-1;

    # create the nodes
    fens = FENodeSetModule.FENodeSet(xyz= xyz);
    # Create the finite elements
    fes = FESetModule.FESetL2(conn=[(1:ncells) (2:ncells+1)]);
    
    return fens,fes;
end
export L2blockx

end

