module FENodeSetModule


using JFinEALE.JFFoundationModule

import Base.count

# Finite element node set class.
type FENodeSet
    # Array of node locations.
    # Array of coordinates, the number of rows corresponds to the number of
    # nodes in the set and the columns corresponds to the space dimensions.
    # The location of node j is given by x(j,:)
    xyz::JFFltMat

    function FENodeSet(xyz::JFFltMat)
        self = new(deepcopy(xyz)); # Need to make a COPY of the input array!
        return self
    end
    
end
export FENodeSet

FENodeSet(;xyz=[]) = FENodeSet(xyz)


# Get the  dimension of the coordinate that defines the location  of the node.
spacedim(self::FENodeSet) =size(self.xyz,2)
export spacedim

function xyz3 (self::FENodeSet)
    # Get the  3-D coordinate that define the location  of the node.
    # Even if the nodes  were specified in  lower dimension (1-D, 2-D)
    # this function returns  a 3-D coordinate  by padding with zeros.
    if (size(self.xyz,2)==1)
        val = [self.xyz zeros(size(self.xyz,1),2)];
    elseif (size(self.xyz,2)==2)
        val = [self.xyz zeros(size(self.xyz,1),1)];
    else
        val = self.xyz;
    end
end
export xyz3

#Get the number of finite element nodes in the node set.
function count(self::FENodeSet)
    return size(self.xyz,1)::JFInt
end
export count


end



