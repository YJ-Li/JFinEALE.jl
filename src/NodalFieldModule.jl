module NodalFieldModule

using JFinEALE.JFFoundationModule

type NodalField{T<:Number}
    name::String
    values::JFMat{T} # Array of degree of freedom parameters,  indexed by node
    dofnums::JFIntMat # Array of degree of freedom numbers, indexed by node
    is_fixed::Matrix{Bool} # Array of Boolean flags, indexed by node
    fixed_values::JFMat{T} # Array of fixed values, indexed by node
    nfreedofs::JFInt # Total number of free degrees of freedom
end
export NodalField

# Constructor of nodal field
function NodalField{T<:Number}(;name="fld",data::JFMat{T}=[],nfens=0,ndofn=1)
    if (length(data)!=0) && (nfens==0)
        values=data
        dofnums=0*similar(values,JFInt,size(values))
        is_fixed=false*similar(values,Bool,size(values))
        fixed_values=zeros(T,size(values))
        nfreedofs=0
    elseif (nfens>0)
        values=zeros(T,nfens,ndofn);
        dofnums=0*similar(values,JFInt,size(values))
        is_fixed=false*similar(values,Bool,size(values))
        fixed_values=zeros(T,size(values))
        nfreedofs=0
    end
    return NodalField(name,values,dofnums,is_fixed,fixed_values,nfreedofs)
end

#Dimension of the nodal parameters (i. e.  how many degrees of freedom per node).
ndofn(self::NodalField) =  size(self.values, 2)
export ndofn

#Number of nodes associated with the field.
nfens(self::NodalField)= size(self.values, 1)
export nfens

#Gather values from the field for the whole system vector.
function gathersysvec{T<:Number} (self::NodalField{T})
    nfens,dim = size(self.values)
    vec=zeros(T,self.nfreedofs,1)
    for i=1:nfens
        for j=1:dim
            en = self.dofnums[i,j]
            if (en > 0)&&(en <= self.nfreedofs)
                vec[en]=self.values[i,j]
            end
        end
    end
    return vec
end
export gathersysvec

# Gather values from the field "as vector".
function gathervaluesasvec!{T<:Number}(self::NodalField,dest::JFVec{T},conn::JFIntVec)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[en]=self.values[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gathervaluesasvec!

# Gather values from the field "as vector".
function gathervaluesasvec!{T<:Number}(self::NodalField,dest::JFVec{T},conn::JFIntMat)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[en]=self.values[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gathervaluesasvec!

# Gather values from the field "as vector".
function gathervaluesasvec!{T<:Number}(self::NodalField,dest::JFMat{T},conn::JFIntMat)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[en]=self.values[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gathervaluesasvec!

# Gather values from the field "as vector".
function gathervaluesasvec!{T<:Number}(self::NodalField,dest::JFMat{T},conn::JFIntVec)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[en]=self.values[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gathervaluesasvec!

# Gather values from the field "as a matrix".
function gathervaluesasmat!{T<:Number}(self::NodalField,dest::JFMat{T},conn::JFIntVec)
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[i,j]=self.values[conn[i],j];
        end
    end
    return self
end
export gathervaluesasmat!

# Gather values from the field "as a matrix".
function gathervaluesasmat!{T<:Number}(self::NodalField,dest::JFMat{T},conn::JFIntMat)
    for i=1:length(conn)
        for j=1:size(self.values,2)
            dest[i,j]=self.values[conn[i],j];
        end
    end
    return self
end
export gathervaluesasmat!

# Gather dofnums from the field "as vector".
function gatherdofnumsasvec!(self::NodalField,dest::JFIntMat,conn::JFIntMat)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.dofnums,2)
            dest[en]=self.dofnums[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gatherdofnumsasvec!

# Gather dofnums from the field "as vector".
function gatherdofnumsasvec!(self::NodalField,dest::JFIntMat,conn::JFIntVec)
    en::JFInt=1;
    for i=1:length(conn)
        for j=1:size(self.dofnums,2)
            dest[en]=self.dofnums[conn[i],j];
            en=en+1;
        end
    end
    return self
end
export gatherdofnumsasvec!

# Number the degrees of freedom.
function numberdofs!(self::NodalField)

    #% The free components in the field are numbered consecutively.
    #% No effort is made to optimize the numbering in any way. If you'd like to
    #% optimize the numbering of the degrees of freedom, use the above form that sets
    #% the permutation of the degrees of freedom, or the permutation of the nodes.

    fixed_dofnum = 0
    nfens,dim = size(self.values)
    self.nfreedofs =0
    for i=1:nfens
        for j=1:dim
            if !self.is_fixed[i,j] # free degree of freedom
                self.nfreedofs = self.nfreedofs + 1
                self.dofnums[i,j] = self.nfreedofs
            else # fixed degree of freedom: no equation
                self.dofnums[i,j] = fixed_dofnum
            end
        end
    end
    return  self
end
export numberdofs!

# Set the EBCs (essential boundary conditions).
function setebc!{T<:Number}(self::NodalField,
                            fenids::JFInt,
                            is_fixed::Bool,
                            comp::JFInt,
                            val::T)
    return setebc!(self, [fenids], [is_fixed], [comp], [val])
end
export setebc!
    
# Set the EBCs (essential boundary conditions).
function setebc!{T<:Number}(self::NodalField,
                            fenids::JFIntVec,
                            is_fixed::Union(BitArray, Array{Bool}),
                            comp::JFIntVec,
                            val::JFVec{T})
    #%
    #% function retobj = set_ebc (self, fenids, is_fixed, comp, val)
    #%
    #% To actually apply the EBCs (that is to copy the fixed values
    #% into the appropriate places in the values array), use apply_ebc().
    #%    Call as:
    #%       retobj = set_ebc (self, fenids, is_fixed, comp, val)
    #%    where
    #%       self           - field
    #%       fenids         - array of N node identifiers
    #%       is_fixed  - an array of N Boolean values which say whether the component
    #%                        is being fixed (value ~=0), or whether
    #%                        it is to be free (value ==0).
    #%       comp           - an array of N components to which to apply the EBCs; comp
    #%                        may be supplied as '[]', in which case the
    #%                        fixed values apply to all components.
    #%       val            - an array of N values of fixed component magnitudes,
    #%                        length(val) == length(fenids)
    #%
    #% Note:  Any call to set_ebc() that would change the current assignment
    #% which degrees of freedom are free and which are fixed invalidates the
    #% current degree-of-freedom numbering. In such a case this method sets
    #%    nfreedofs = 0; and
    #%    dofnums=dofnums*0+sysmat_assembler.invalid_dofnum;
    #%
    Maxl =max(length(fenids), length(is_fixed), length(comp), length(val))
    if (length(fenids)==1)&&(Maxl>1)
        fenids=repmat(fenids,Maxl,1);
    end
    if (length(is_fixed)==1)&&(Maxl>1)
        is_fixed=repmat(is_fixed,Maxl,1);
    end
    if (length(comp)==1)&&(Maxl>1)
        comp=repmat(comp,Maxl,1);
    end
    if (length(val)==1)&&(Maxl>1)
        val=repmat(val,Maxl,1);
    end
    nfens,dim = size(self.values)
    for i=1:length(fenids)
        k = fenids[i];
        if (isempty(comp))
            m = 1:dim;
        else
            m = comp[i];
        end
        self.is_fixed[k,m] = is_fixed[i];
        if self.is_fixed[k,m]
            self.fixed_values[k,m] = val[i];
        end
    end
    #    % Any call to set_ebc() which changes which degrees of
    #    %     freedom are free and which are fixed
    #    %     invalidates the current equation numbering.
    self.nfreedofs = 0
    self.dofnums=0*self.dofnums
    return  self
end
export setebc!

# Apply EBCs (essential boundary conditions).
function applyebc!(self::NodalField)
    #%
    #% function retobj = apply_ebc (self)
    #%
    #% This function applies existing boundary conditions.  In other
    #% words, it sets value=fixed_values for all parameters and components
    #% of parameters that are being fixed as essential boundary
    #% conditions (EBCs).
    #%
    #% To set the prescribed values and to mark components as being prescribed,
    #% use set_ebc().
    #%
    #% Don't forget to assign the field to itself, because this function
    #% changes the internal state of self.
    #%    Call as:
    #%      f = apply_ebc(f);
    #%
    nfens,dim = size(self.values);
    for i=1:nfens
        for j=1:dim
            if self.is_fixed[i,j]
                self.values[i,j] = self.fixed_values[i,j];
            end
        end
    end
    return  self
end
export applyebc!

# Scatter values to the field from a system vector.
function scattersysvec!{T<:Number}(self::NodalField, vec::JFVec{T})
    #%
    #% function retobj = scatter_sysvec (self, vec)
    #%
    #%   Call as:
    #%     f=scatter_sysvec(f, vec)
    #%   where
    #%     f=field (note that we have to assign to f: f changes inside this function)
    #%     vec=system vector (an array of doubles), indexed by free degree of freedom numbers
    #%
    nfens,dim = size(self.values);
    for i=1:nfens
        for j=1:dim
            dn = self.dofnums[i,j];
            if (dn > 0) && (dn <= self.nfreedofs)
                self.values[i,j] = vec[dn];
            end
        end
    end
    return  self
end
export scattersysvec!

function  clone{T<:Number}(self::NodalField{T}, name:: String="")
# % Make a new nodal field by cloning another field
# %
# % function retobj = clone(fld,name)
# %
# %   Call as:
# %     clonedf = clone(f, name)
# %   where
# %     f=field to clone
# %     name=name of the new field
# %
    retobj =  deepcopy(self);
    retobj.name = name;
    return retobj;
end
export clone


end

