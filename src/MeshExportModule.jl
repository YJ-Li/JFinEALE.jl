module MeshExportModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.FENodeSetModule

const L2=3
const T3=5
const Q4=9
const T4=10
const H8=12
const Q8=23

# Export mesh to a VTK 1.0 file as an unstructured grid.
function vtkexportmesh{T<:FESet} (theFile::String, fens::FENodeSet, fes::T; opts...)
    if typeof(fes)==FESetL2
        Cell_type=L2
    elseif typeof(fes)==FESetT3
        Cell_type=T3
    elseif typeof(fes)==FESetQ4
        Cell_type=Q4
    elseif typeof(fes)==FESetT4
        Cell_type=T4
    elseif typeof(fes)==FESetH8
        Cell_type=H8
    elseif typeof(fes)==FESetQ8
        Cell_type=Q8
    else
        error("Cannot handle $(typeof(fes))")
    end
    return vtkexportmesh (theFile,fes.conn,fens.xyz, Cell_type; opts...)
end

# Export mesh to a VTK 1.0 file as an unstructured grid.
function vtkexportmesh (theFile::String,Connectivity,Points, Cell_type;
                        vectors=nothing,vectors_name ="vectors",
                        scalars=nothing, scalars_name ="scalars",binary = false)
# theFile= File name as a string,
    # Connectivity= connectivity array, one row per cell,
    # Points= Coordinate array, one row per point,
    # Cell_type= Cell type code: L2=3, T3=5, Q4=9, T4=10, H8=12
    # options = structure,  with optional attributes:
    # options.binary= Boolean flag, should the file be written as binary?
    # Default is false (file is ASCII).
    # options.scalars =Array of per-node data, same number of rows as  array Points.
    # options.scalars_name = String (no spaces), name for the scalars array.
    # options.vectors =Array of per-node data, same number of rows as  array Points, three columns.
    # options.vectors_name = String (no spaces), name for the vectors array.

    X= Points
    if size(Points, 2)<3
        X=   zeros(size (Points,1),3)
        X [:,1: size(Points,2)] = Points
    end
    if (Cell_type==L2) && (size(Connectivity,2)!=2)
        error("Wrong number of connected nodes")
    end
    if (Cell_type==T3) && (size(Connectivity,2)!=3)
        error("Wrong number of connected nodes")
    end
    if (Cell_type==Q4) && (size(Connectivity,2)!=4)
        error("Wrong number of connected nodes")
    end
    if (Cell_type==T4) && (size(Connectivity,2)!=4)
        error("Wrong number of connected nodes")
    end
    if (Cell_type==H8) && (size(Connectivity,2)!=8)
        error("Wrong number of connected nodes")
    end

    #    if (~strcmp(ext,"vtk"))
    #        theFile = [theFile ".vtk"];
    #    end

    fid=open(theFile,"w");
    if (fid==-1)
        error (["Could not open " * theFile])
        return nothing
    end
    print(fid,"# vtk DataFile Version 1.0\n");
    print(fid,"JFinEALE Export\n");
    if (!binary)
        print(fid,"ASCII\n");
    else
        print(fid,"BINARY\n");
    end
    print(fid,"\n");
    print(fid,"DATASET UNSTRUCTURED_GRID\n");
    print(fid,"POINTS ", size(X,1), " double\n");
    if (!binary)
        for i= 1:size(X, 1)
            for j= 1:size(X,2)-1
                print(fid,X[i,j]," ");
            end
            print(fid,X[i,end],"\n");
        end
    else
        fwrite(fid,cast(X,"double"),"double","n");
    end
    print(fid,"\n");
    print(fid,"\n");

    print(fid,"CELLS ",size(Connectivity,1)," ",(size(Connectivity,1)*(size(Connectivity,2)+1)),"\n");
    if (!binary)
        for i= 1:size(Connectivity, 1)
            print(fid,size(Connectivity,2)," ");
            for j= 1:size(Connectivity,2)-1
                print(fid,Connectivity[i,j]-1," ");
            end
            print(fid,Connectivity[i,end]-1,"\n");
        end
    else
        fwrite(fid,cast([zeros(size(f,1),1)+size(f,2),f-1],"int32"),"int32","n");
    end
    print(fid,"\n");
    print(fid,"\n");
    print(fid,"CELL_TYPES ",size(Connectivity,1),"\n");
    if (!binary)
        for i= 1:size(Connectivity,1)
            print(fid,Cell_type,"\n");
        end
    else
        fwrite(fid,cast(ctype,"int32"),"int32","n");
    end

    if scalars!=nothing
        print(fid,"POINT_DATA ",length(scalars),"\n");
        print(fid,"SCALARS ",scalars_name," double\n");
        print(fid,"LOOKUP_TABLE default\n");
        if (!binary)
            for j= 1:length(scalars)
                print(fid,scalars[j],"\n");
            end
        else
            fwrite(fid,cast(scalars,"double"),"double","n");
        end
    end

    if vectors!=nothing
        print(fid,"POINT_DATA ",size(vectors,1),"\n");
        print(fid,"VECTORS ",vectors_name," double\n");
        #print(fid,"LOOKUP_TABLE default\n");
        if size(vectors, 2)<3
            X=   zeros(size (vectors,1),3)
            X [:,1: size(vectors,2)] = vectors
        end
        if (!binary)
            for j= 1:size(X,1)
                k=1;
                print(fid,X[j,k]);
                for k=2:size(X,2)
                    print(fid," ",X[j,k]);
                end
                print(fid,"\n");
            end
        else
            fwrite(fid,cast(X,"double"),"double","n");
        end
    end

    print(fid,"\n");
    fid=close(fid);
    return true
end
export vtkexportmesh

function finealemesh(fens,fes,meshmfile)
    fid=open(meshmfile * ".m","w");
    if (fid==-1)
        error("Could not open " * meshmfile)
        return nothing
    end
    print(fid,"function [fens,fes]= $meshmfile()\n");
    print(fid,"fens=fenode_set(struct('xyz',[\n");
    for i= 1:size(fens.xyz, 1)
        for j= 1:size(fens.xyz,2)-1
            print(fid,fens.xyz[i,j]," ");
        end
        print(fid,fens.xyz[i,end],"\n");
    end
    print(fid,"]));\n");
    # write out the finite element set
    s=string(typeof(fes));
    m=match(r"Module", s)
    if (m!=nothing)
        s=s[m.offset+6:end]
    end
    m=match(r"FESet", s)
    T=s[m.offset+5:end]         # type of the finite element
    print(fid,"fes =fe_set_$(T)(struct ('conn',[\n");
    for i= 1:size(fes.conn, 1)
        for j= 1:size(fes.conn,2)-1
            print(fid,fes.conn[i,j]," ");
        end
        print(fid,fes.conn[i,end],"\n");
    end
    print(fid,"]));\n");
    print(fid,"end\n");
    fid=close(fid);
end
export finealemesh

function graphcsv(xdata,ydata;xlabel="x",ylabel="y",file= "agraph")
    fid=open(file * ".csv","w");
    if (fid==-1)
        error("Could not open " * file)
        return nothing
    end
    print(fid,"$(xlabel),$(ylabel)\n");
    for i= 1:length(xdata)
        print(fid,"$(xdata[i]),$(ydata[i])\n");
    end
    fid=close(fid);
end
export graphcsv


end


