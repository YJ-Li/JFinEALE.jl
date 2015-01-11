module MeshTransformModule


using JFFoundationModule
using FENodeSetModule
using FESetModule

function mirrormesh{T<:FESet}(fens::FENodeSet, fes::T, Normal::JFFltVec, Point::JFFltVec; args...)
    # % Mirror a 2-D mesh in a plane given by its normal and one point. 
    # %
    # % function [fens1,fes1] = mirror_mesh(fens, fes, Normal, Point, renumb)
    # %
    # % fens, gcells= mesh, 
    # % Normal, Point= 2-D arrays
    # %  
    # % Warning: The code to relies on the numbering of the cells: to reverse
    # % the orientation of the mirrored cells, the connectivity is listed in
    # % reverse order.   If the mirrored cells do not follow this rule (for instance 
    # % hexahedra for quadrilaterals), their areas will
    # % come out negative. In such a case the renumbering function 
    # % of the connectivity needs to be supplied..
    # %
    # % For instance: H8 elements require  the renumbering function to be supplied as
    # % [fens1,gcells1] = mirror_mesh(fens, gcells,...
    # %           [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    # %

    # Treat optional arguments
    renumb(conn)=conn[end:-1:1]; # Simply switch the order of nodes.  Works for simplexes...
    for arg in args
        sy, val = arg
        if sy==:renumb
            renumb=val
        end
    end
    # Make sure we're doing with a unit normal
    Normal = Normal/norm (Normal);

    fens1 =deepcopy(fens); # the mirrored mesh nodes
    for i=1:count(fens1)
        a=fens1.xyz[i,:]
        d=(Normal'*vec(a-reshape(Point,1,length(Point))));
        fens1.xyz[i,:] =a-2*d*reshape(Normal,1,size(fens1.xyz,2));
    end
    # Reconnect the cells
    fes1=deepcopy(fes);
    for i=1:size(fes1.conn,1)
        fes1.conn[i,:]=renumb(fes1.conn[i,:]);
    end
    return fens1,fes1 
end
export mirrormesh

end


