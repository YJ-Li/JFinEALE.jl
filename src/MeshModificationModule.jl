module MeshModificationModule

using JFinEALE.JFFoundationModule
using JFinEALE.FENodeSetModule
using JFinEALE.FESetModule
using Base.Sort
using Base.Order

function meshboundary{T<:FESet}(fes::T)
    # Extract the boundary finite elements from a mesh.
    #
    # function bdry_fes = mesh_boundary(fes, options)
    #
    # Extract the finite elements of manifold dimension (n-1) from the
    # supplied list of finite elements of manifold dimension (n).
    #    options = struct with any attributes that should be passed to the
    #    constructor of the boundary finite elements
    #

    # Form all hyperfaces, non-duplicates are boundary cells
    hypf= FESetModule.boundaryconn(fes);    # get the connectivity of the boundary elements
    bdryconn =myunique2(hypf);
    make = FESetModule.boundaryfe(fes);     # get the function that can make a boundary element

    return make(conn =bdryconn);
end
export meshboundary

function mysortrows(A::JFIntMat)
    # Sort the rows of A by sorting each column from back to front.

    m,n = size(A);
    
    indx =  zeros(JFInt,m); sindx = zeros(JFInt,m)
    for i=1:m
        indx[i]=i
    end
    nindx =  zeros(JFInt,m); 
    col = zeros(JFInt,m)
    for c = n:-1:1
        for i=1:m
            col[i]=A[indx[i],c]
        end
        #Sorting a column vector is much faster than sorting a column matrix
        sindx=sortperm(col,alg=QuickSort);
        #sortperm!(sindx,col,alg=QuickSort); # available for 0.4, slightly faster
        #indx=indx[sindx] # saving allocations by using the below loops
        for i=1:m
            nindx[i]=indx[sindx[i]]
        end
        for i=1:m
            indx[i]=nindx[i]
        end
    end
        
    return A[indx,:]
end

function mysortdim2!(A::JFIntMat)
    # Sort each row  of A in ascending order.

    m,n = size(A);
    r = zeros(JFInt,n)
   @inbounds for k = 1:m
        for i=1:n
            r[i]=A[k,i]
        end
        sort!(r);
        for i=1:n
            A[k,i]=r[i]
        end
    end
    return A
end

function  myunique2(A::JFIntMat) # speeded up; now the bottleneck is mysortrows
    #println("size(A)=$(size(A))")
    maxA=maximum(A[:])::JFInt
    sA=deepcopy(A)
    @time sA=mysortdim2!(sA)::JFIntMat;#this is fast
    #@time sA=sort(A,2,alg=QuickSort)::JFIntMat;#this is slow
    sA= [sA (1:size(A,1))+maxA]::JFIntMat
    @time sA =mysortrows(sA); # this now takes the majority of time, but much less than the function below
    #@time sA  = sortrows(sA,alg=QuickSort);;#this is slow
    rix=sA[:,end]; rix=rix[:]-maxA;
    sA=sA[:,1:end-1];
    d=falses(size(sA,1)-1)
    for k=1:length(d)
        for m=1:size(sA,2)
            if sA[k,m]!=sA[k+1,m]
                d[k]=true;
                break;
            end
        end
    end
    #d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
    ad=zeros(JFInt,size(d,1)+1)
    ad[1]=1;
    for k=2:length(ad)
        for m=1:size(d,2)
            if d[k-1,m]!=0
                ad[k]=1;
                break;
            end
        end
    end
    #ad=map((x) -> (x?1:0),[true; any(d,2)]);
    iu=trues(length(ad))
    for k=1:(length(ad)-1)
        ad[k]=ad[k]+ad[k+1]
        iu[k]=(ad[k]>1)
    end
    ad[end]=ad[end]+1;
    iu[end]=(ad[end]>1)
    #iu =map((x) -> (x>1? true: false),(ad + [ad[2:end];1]));
    Out =A[rix[iu],:];
    return Out
end

### This code is correct, but very slow.
function  myunique1(A::JFIntMat)
    maxA=maximum(A[:])
    sA=sort(A,2);# most time spent here
    sA= [sA (1:size(A,1))+maxA]
    sA  = sortrows(sA);;#and here
    rix=sA[:,end]; rix=rix[:]-maxA;
    sA=sA[:,1:end-1];
    d=(sA[1:end-1,:].!=sA[2:end,:]); # element-wise comparison!
    ad=map((x) -> (x?1:0),[true; any(d,2)]);
    iu =map((x) -> (x>1? true: false),(ad + [ad[2:end];1]));
    Out =A[rix[iu[:]],:];
    return Out
end


function fusenodes(fens1::FENodeSet, fens2::FENodeSet, tolerance:: JFFlt)
    # % Fuse together nodes from two node sets.
    # %
    # % function [fens,new_indexes_of_fens1_nodes] = fuse_nodes(fens1, fens2, tolerance)
    # %
    # % Fuse two node sets. If necessary, by gluing together nodes located within tolerance of each other.
    # % The two node sets, fens1 and fens2,  are fused together by
    # % merging the nodes that fall within a box of size "tolerance".
    # % The merged node set, fens, and the new  indexes of the nodes
    # % in the set fens1 are returned.
    # %
    # % The set fens2 will be included unchanged, in the same order,
    # % in the node set fens.
    # %
    # % The indexes of the node set fens1 will have changed.
    # %
    # % Example:
    # % After the call to this function we have
    # % k=new_indexes_of_fens1_nodes(j) is the node in the node set fens which
    # % used to be node j in node set fens1.
    # % The finite element set connectivity that used to refer to fens1
    # % needs to be updated to refer to the same nodes in  the set fens as
    # %     fes = update_conn(fes ,new_indexes_of_fens1_nodes);
    # %
    # % See also: merge_nodes, update_conn
    # %
    xyz1::JFFltMat = copy(fens1.xyz);
    id1::JFIntVec =zeros(JFInt,size(xyz1,1));
    for rx=1:length(id1)
        id1[rx]=rx;
    end
    dim =size(xyz1,2);
    xyz2::JFFltMat = copy(fens2.xyz);
    id2::JFIntVec =zeros(JFInt,size(xyz2,1));
    for rx=1:length(id2)
        id2[rx]=rx;
    end
    n1::JFFlt= 0.0
    # % Mark nodes from the first array that are duplicated in the second
    if (tolerance>0.0) # should we attempt to merge nodes?
        for i=1:size(xyz1,1)
            for rx=1:size(xyz2,1)
                n1= 0.0
                for cx=1:size(xyz2,2)
                    n1=n1+abs(xyz2[rx,cx]-xyz1[i,cx]);
                end
                if (n1<tolerance)
                    id1[i] =-rx; break;
                end
            end
        end
    end
    # % Generate  fused arrays of the nodes
    xyzm = zeros(JFFlt,size(xyz1,1)+size(xyz2,1),dim);
    for rx=1:size(xyz2,1)
        xyzm[rx,:]=xyz2[rx,:];
    end
    idm = zeros(JFInt,size(xyz1,1)+size(xyz2,1));
    for rx=1:size(xyz2,1)
        idm[rx]=rx;
    end
    mid=size(xyz2,1)+1;
    for i=1:size(xyz1,1) # and then we pick only non-duplicated fens1
        if id1[i]>0
            id1[i]=mid;
            idm[mid]=mid;
            for cx=1:size(xyz1,2)
                xyzm[mid,cx]=xyz1[i,cx];
            end
            mid=mid+1;
        else
            id1[i]=id2[-id1[i]];
        end
    end
    nfens =mid-1;
    xyzm =xyzm[1:nfens,:];

    # % Create the fused Node set
    fens =FENodeSetModule.FENodeSet(xyz=xyzm);
    # % The Node set 1 numbering will change
    new_indexes_of_fens1_nodes=id1[:];
    # % The node set 2 numbering stays the same
    return fens, new_indexes_of_fens1_nodes
end
export fusenodes


function compactfens(fens::FENodeSetModule.FENodeSet, connected::JFIntVec)
    # % Compact the finite element node set by deleting unconnected nodes.
    #
    # % fens = array of finite element nodes
    # % connected = The array element connected(j) is either 0 (when j is an unconnected
    # %    node), or a positive number (when node j is connected to other nodes by
    # %    at least one finite element)
    # %
    # % Output:
    # % fens = new set of finite element nodes
    # % new_numbering= array which tells where in the new fens array the
    # %      connected nodes are (or 0 when the node was unconnected). For instance,
    # %      node 5 was connected, and in the new array it is the third node: then
    # %      new_numbering(5) is 3.
    # %
    # % Examples:
    # %
    # % Let us say there are nodes not connected to any finite element that you
    # % would like to remove from the mesh: here is how that would be
    # % accomplished.
    # %
    # % connected = find_unconn_fens(fens, fes);
    # % [fens, new_numbering] =compact_fens(fens, connected);
    # % fes = renumber_fe_conn(fes, new_numbering);
    # %
    # % Finally, check that the mesh is valid:
    # % validate_mesh(fens, fes);
    # %

    new_numbering=zeros(JFInt,count(fens),1);
    xyz =fens.xyz;
    nxyz =xyz;
    id=1;
    for i=1:length(connected)
        if (connected[i]>=0)
            new_numbering[i]=id;
            nxyz[id,:] =xyz [i,:];
            id=id+1;
        end
    end
    new_numbering=new_numbering[1:id-1];
    fens=FENodeSetModule.FENodeSet (xyz=nxyz[1:id-1,:]);
    return fens, new_numbering[1:id-1]
end
export compactfens

function mergemeshes{T<:FESet}(fens1::FENodeSet, fes1::T, fens2::FENodeSet, fes2::T, tolerance::JFFlt)
    # % Merge together two meshes.
    # %
    # % function [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens2,
    # %                     fes2, tolerance)
    # %
    # % Merge two meshes together by gluing together nodes within tolerance. The
    # % two meshes, fens1, fes1, and fens2, fes2, are glued together by merging
    # % the nodes that fall within a box of size "tolerance". If tolerance is set
    # % to zero, no merging of nodes is performed; the two meshes are simply
    # % concatenated together.
    # %
    # % The merged node set, fens, and the two arrays of finite elements with
    # % renumbered  connectivities are returned.
    # %
    # % Important notes: On entry into this function the connectivity of fes1
    # % point into fens1 and the connectivity of fes2 point into fens2. After
    # % this function returns the connectivity of both fes1 and fes2 point into
    # % fens. The order of the nodes of the node set fens1 in the resulting set
    # % fens will have changed, whereas the order of the nodes of the node set
    # % fens2 is are guaranteed to be the same. Therefore, the connectivity of
    # % fes2 will in fact remain the same.
    # %
    # %
    # % See also: fusenodes, updateconn
    # %

    # % Fuse the nodes
    #show(code_typed( fusenodes,(FENodeSet, FENodeSet, JFFlt)))
    fens,new_indexes_of_fens1_nodes = fusenodes(fens1, fens2, tolerance);
    # % Renumber the finite elements
    newfes1= deepcopy(fes1)
    updateconn!(newfes1,new_indexes_of_fens1_nodes);
    # % Note that now the connectivity of both fes1 and fes2 point into
    # % fens.
    return fens,newfes1,fes2
end
export mergemeshes

function mergenmeshes!(fensa, fesa, tolerance::JFFlt)
#     % Merge several meshes together.
# %
# % function [fens,fesa] = merge_n_meshes(fensa, fesa, tolerance)
# %
# % Merge several meshes together either by simple concatenation of nodes or by
# % gluing together nodes within tolerance.
# %
# % Inputs:
# % fensa= cell array of node sets, one for each mesh;
# % fesa= cell array of finite element sets, one for each mesh;
# % tolerance= Geometric tolerance, maybe supplied as zero (>=0).
# %
# % The meshes are glued together by
# % merging the nodes that fall within a box of size "tolerance". If tolerance
# % is set to zero, no merging of nodes is performed; the nodes from the meshes are
# % simply concatenated together.
# %
# % The merged node set, fens, and the cell array of finite element sets with
# % renumbered  connectivities are returned.
# %
# % Outputs:
# % fens= merged node set,
# % fesa= cell array of finite element sets updated to use the merged node set.
# %
# %
# % See also: merge_meshes


    if (length(fensa))!=(length(fesa))
        error("(length(fensa))!=(length(fesa))");
    end
    if (length(fensa))==1
        fens=fensa[1];
        return fens,fesa        # There is nothing to be done: this is a single mesh
    end
    fens=fensa[1];
    for j=2:length(fesa)
        fens,new_indexes_of_fens1_nodes = fusenodes(fensa[j], fens, tolerance);
        FESetModule.updateconn!(fesa[j],new_indexes_of_fens1_nodes);
    end
    return fens,fesa
end
export mergenmeshes!

function mergenodes(fens::FENodeSetModule.FENodeSet, fes::FESetModule.FESet, tolerance::JFFlt)
    # % Merge together  nodes of a single node set.
    # %
    # % function [fens,fes] = merge_nodes(fens, fes, tolerance)
    # %
    # % Merge by gluing together nodes from a single node set located within
    # % tolerance of each other. The nodes are glued together by merging the
    # % nodes that fall within a box of size "tolerance". The merged node
    # % set, fens, and the finite element set with renumbered  connectivities
    # % are returned.
    # %
    # % See also: fuse_nodes

    xyz1 = fens.xyz;
    dim =size(xyz1,2);
    id1 = (1:FENodeSetModule.count(fens))';
    c1=ones(size(xyz1,1),1);
    xyzd= zeros(size(xyz1));
    d= zeros(size(xyz1,1));
    m= trues(size(xyz1,1));
    # Mark nodes from the array that are duplicated
    for i=1:FENodeSetModule.count(fens)
        if (id1[i]>0) # This node has not yet been marked for merging
            XYZ =xyz1[i,:];
            xyzd[:,:]=abs(xyz1-c1*XYZ); #find the distances along  coordinate directions
            d=sum(xyzd,2);
            map!((x)->x<tolerance, m, d);
            jx=find(m);
            if (!isempty(jx))
                minn=minimum(jx);
                id1[jx] =-minn;
                id1[minn] =minn;
            end
        end
    end
    # Generate  merged arrays of the nodes
    xyzm = zeros(JFFlt,FENodeSetModule.count(fens),dim);
    mid=1;
    for i=1:FENodeSetModule.count(fens) # and then we pick only non-duplicated fens1
        if id1[i]>0
            id1[i]=mid;
            xyzm[mid,:]=xyz1[i,:];
            mid=mid+1;
        else
            id1[i]=id1[-id1[i]];
        end
    end
    nfens =mid-1;
    xyzm =xyzm[1:nfens,:];
    # Renumber the cells
    conns=fes.conn;
    for i=1:FESetModule.count(fes)
        conn=conns[i,:];
        conns[i,:]=id1[conn];
    end
    fes.conn=conns;

    fens=FENodeSetModule.FENodeSet(xyz=xyzm[1:nfens,:]);

    return fens,fes
end
export mergenodes

function renumberfeconn(fes::FESetModule.FESet, new_numbering::JFIntVec)
    # % Renumber the nodes in the connectivity of the finite elements based on a new
    # % numbering for the nodes.
    # %
    # % function fes = renumber_fe_conn(fes, new_numbering)
    # %
    # % fes =finite element set
    # % new_numbering = new serial numbers for the nodes.  The connectivity
    # %           should be changed as conn(j) --> new_numbering(conn(j))
    # %
    # % Output:
    # % fes = new Finite element set
    # %
    # % Let us say there are nodes not connected to any finite element that you
    # % would like to remove from the mesh: here is how that would be
    # % accomplished.
    # %
    # % connected = find_unconn_fens(fens, fes);
    # % [fens, new_numbering] =compact_fens(fens, connected);
    # % fes = renumber_fe_conn(fes, new_numbering);
    # %
    # % Finally, check that the mesh is valid:
    # % validate_mesh(fens, fes);
    # %
    conn =fes.conn;
    for i=1:size(conn,1)
        c=conn[i,:];
        conn[i,:] =new_numbering[c]';
    end
    fes.conn=conn;
    return fes
end
export renumberfeconn

end



