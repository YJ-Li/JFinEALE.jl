module MeshTriangleModule

using JFFoundationModule
using FESetModule
using FENodeSetModule

# T3 Mesh of a rectangle.
function T3block (Length::JFFlt,Width::JFFlt,nL::JFInt,nW::JFInt)
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

    nnodes=(nL+1)*(nW+1);
    ncells=2*(nL)*(nW);
    xs=zeros(JFFlt,nnodes,2);
    conns=zeros(JFInt,ncells,3);
    f=1;
    for j=1:(nW+1)
        for i=1:(nL+1)
            xs[f,:]=[(i-1)*Length/nL (j-1)*Width/nW];
            f=f+1;
        end
    end
    fens=FENodeSetModule.FENodeSet(xyz=xs);

    function node_numbers1(i,j,nL,nW)
        f=(j-1)*(nL+1)+i;
        nn=[f (f+1) f+(nL+1)];
    end
    function node_numbers2(i,j,nL,nW)
        f=(j-1)*(nL+1)+i;
        nn=[(f+1)  f+(nL+1)+1 f+(nL+1)];
    end

    gc=1;
    for i=1:nL
        for j=1:nW
            nn=node_numbers1(i,j,nL,nW);
            conns[gc,:]=nn;
            gc=gc+1;
            nn=node_numbers2(i,j,nL,nW);
            conns[gc,:]=nn;
            gc=gc+1;
        end
    end
    fes = FESetModule.FESetT3(conn=conns);
    
    return fens,fes
end
export T3block


end

