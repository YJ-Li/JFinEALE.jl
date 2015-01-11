module AssemblyModule

using JFFoundationModule

const  inv_dofnum=0;            # invalid degree of freedom number -- no equation

abstract SysmatAssemblerBase;
export SysmatAssemblerBase

type SysmatAssemblerSparse{T<:Number} <: SysmatAssemblerBase
  # Type for assembling of a sparse global matrix from elementwise matrices.
  buffer_length::JFInt;
  matbuffer::JFVec{T};
  rowbuffer::JFIntVec;
  colbuffer::JFIntVec;
  buffer_pointer::JFInt;
  ndofs_row::JFInt; ndofs_col::JFInt;
end
export SysmatAssemblerSparse

# AssemblyModule.SysmatAssemblerSparse()
function SysmatAssemblerSparse{T<:Number}(zero::T=0.0)
    return SysmatAssemblerSparse{T}(0,[zero],[0],[0],0,0,0)
end

"""
Start the assembly of a global matrix.
"""
function startassembly!{T<:Number}(self::SysmatAssemblerSparse{T},
                        elem_mat_nrows::JFInt, elem_mat_ncols::JFInt, elem_mat_nmatrices::JFInt,
                        ndofs_row::JFInt, ndofs_col::JFInt)
  # Start the assembly of a global matrix.
  # The method makes buffers for matrix assembly. It must be called before
  # the first call to the method assemble.
  # elem_mat_nrows= number of rows in typical element matrix,
  # elem_mat_ncols= number of columns in a typical element matrix,
  # elem_mat_nmatrices= number of element matrices,
  # ndofs_row= Total number of equations in the row direction,
  # ndofs_col= Total number of equations in the column direction.
  self.buffer_length =elem_mat_nmatrices*elem_mat_nrows*elem_mat_ncols;
  self.rowbuffer = zeros (JFInt,self.buffer_length);
  self.colbuffer = zeros (JFInt,self.buffer_length);
  self.matbuffer = zeros (T,self.buffer_length);
  self.buffer_pointer = 1;
  self.ndofs_row =ndofs_row;
  self.ndofs_col =ndofs_col;
  return self
end
export startassembly!

function assemble!{T<:Number}(self::SysmatAssemblerSparse{T}, mat::JFMat{T},
                              dofnums_row::JFIntMat, dofnums_col::JFIntMat)
  # Assembly of a rectangular matrix.
  # The method assembles a rectangular matrix using the two vectors of
  # equation numbers for the rows and columns.
  nrows=length(dofnums_row); ncolumns=length(dofnums_col);
  p = self.buffer_pointer
  for j=1:ncolumns
    for i=1:nrows
      self.matbuffer[p] = mat[i,j] # serialized matrix
      self.rowbuffer[p] = dofnums_row[i];
      self.colbuffer[p] = dofnums_col[j];
      p=p+1
    end
  end
  self.buffer_pointer=p;
  return self
end
export assemble!


function makematrix!(self::SysmatAssemblerSparse)
  # Make a sparse matrix.
  # The method makes a sparse matrix from the assembly buffers.
  for j=1:self.buffer_pointer-1
    if self.rowbuffer[j] == inv_dofnum
      self.rowbuffer[j]=self.ndofs_row+1;
    end
    if self.colbuffer[j] == inv_dofnum
      self.colbuffer[j]=self.ndofs_col+1;
    end
  end
  S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
             self.colbuffer[1:self.buffer_pointer-1],
             self.matbuffer[1:self.buffer_pointer-1],
             self.ndofs_row+1, self.ndofs_col+1);
  self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
  return S[1:end-1,1:end-1]
end
export makematrix!





type SysmatAssemblerSparseSymm{T<:Number} <: SysmatAssemblerBase
  # Type for assembling of a sparse global matrix from elementwise matrices.
  buffer_length:: JFInt;
  matbuffer::JFVec{T};
  rowbuffer::JFIntVec;
  colbuffer::JFIntVec;
  buffer_pointer:: JFInt;
  ndofs_row:: JFInt; ndofs_col:: JFInt;
end
export SysmatAssemblerSparseSymm

# AssemblyModule.SysmatAssemblerSparse()
function SysmatAssemblerSparseSymm{T<:Number}(zero::T=0.0)
    return SysmatAssemblerSparseSymm{T}(0,[zero],[0],[0],0,0,0)
end

"""
Start the assembly of a global matrix.
"""
function startassembly!{T<:Number}(self::SysmatAssemblerSparseSymm{T},
                        elem_mat_nrows::JFInt, elem_mat_ncols::JFInt, elem_mat_nmatrices::JFInt,
                        ndofs_row::JFInt, ndofs_col::JFInt)
  # Start the assembly of a global matrix.
  # The method makes buffers for matrix assembly. It must be called before
  # the first call to the method assemble.
  # elem_mat_nrows= number of rows in typical element matrix,
  # elem_mat_ncols= number of columns in a typical element matrix,
  # elem_mat_nmatrices= number of element matrices,
  # ndofs_row= Total number of equations in the row direction,
  # ndofs_col= Total number of equations in the column direction.
  self.buffer_length =elem_mat_nmatrices*elem_mat_nrows*elem_mat_ncols;
  self.rowbuffer = zeros (JFInt,self.buffer_length);
  self.colbuffer = zeros (JFInt,self.buffer_length);
  self.matbuffer = zeros (T,self.buffer_length);
  self.buffer_pointer = 1;
  self.ndofs_row =ndofs_row;
  self.ndofs_col =ndofs_col;
  return self
end
export startassembly!

function assemble!{T<:Number}(self::SysmatAssemblerSparseSymm{T}, mat::JFMat{T},
                              dofnums_row::JFIntMat, dofnums_col::JFIntMat)
  # Assembly of a rectangular matrix.
  # The method assembles a rectangular matrix using the two vectors of
  # equation numbers for the rows and columns.
  nrows=length(dofnums_row); ncolumns=length(dofnums_col);
  p = self.buffer_pointer
  for j=1:ncolumns
    for i=j:nrows
      self.matbuffer[p] = mat[i,j] # serialized matrix
      self.rowbuffer[p] = dofnums_row[i];
      self.colbuffer[p] = dofnums_col[j];
      p=p+1
    end
  end
  self.buffer_pointer=p;
  return self
end
export assemble!


function makematrix!(self::SysmatAssemblerSparseSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    for j=1:self.buffer_pointer-1
        if self.rowbuffer[j] == inv_dofnum
            self.rowbuffer[j]=self.ndofs_row+1;
        end
        if self.colbuffer[j] == inv_dofnum
            self.colbuffer[j]=self.ndofs_col+1;
        end
    end
    S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs_row+1, self.ndofs_col+1);
    S = S+S';    # construct the other triangle
    for j=1:size(S,1)
        S[j,j]=S[j,j]/2.0;      # the diagonal is there twice; fix it;
    end
    self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
    return S[1:end-1,1:end-1]
end
export makematrix!



abstract SysvecAssemblerBase;
export SysvecAssemblerBase

type SysvecAssembler{T<:Number} <: SysvecAssemblerBase
  # % The class sysvec_assembler is for assembling of a system
  # % column vector from elementwise vectors.
  # 	% Fixed degrees of freedom numbers are given this value:
  #      % it indicates that this is not a valid free  degree of freedom number.
  F_buffer::JFVec{T};
end
export SysvecAssembler

function SysvecAssembler{T<:Number}(zero::T=0.0)
    return SysvecAssembler([zero])
end


function startassembly!{T<:Number}(self::SysvecAssembler{T}, ndofs_row::JFInt)
  # Start assembly.
  # 
  # The method makes the buffer for the vector assembly. It must be called before
  #  the first call to the method assemble.
  #  ndofs_row= Total number of degrees of freedom.
  self.F_buffer= zeros(T,ndofs_row);
end
export startassembly!

function assemble!{T<:Number}(self::SysvecAssembler{T}, vec::JFMat{T}, dofnums::JFIntMat)
  #   % Assembly of elementwise vector.
  #   %
  #   % function assemble(self, vec, dofnums)
  #     %
  #     % The method assembles a column element vector using the vector of
  #     % equation numbers for the rows.
  for i = 1:length(dofnums)
    gi = dofnums[i];
    if (gi != inv_dofnum)
      self.F_buffer[gi] = self.F_buffer[gi] + vec[i];
    end
  end
end
export assemble!

"""
Docu
"""
function makevector! (self::SysvecAssembler)
  # Make the global vector.
  # The method makes a vector from the assembly buffers.
  return deepcopy(self.F_buffer);
end
export makevector!


end
