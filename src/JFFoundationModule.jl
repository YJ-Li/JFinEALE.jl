module JFFoundationModule

# Documentation relies on this
#VERSION < v"0.4-" && using Docile

import Base.Complex

typealias JFInt Int
typealias JFFlt Float64
typealias JFCplxFlt Complex{Float64}
typealias JFFltVec Vector{JFFlt}
typealias JFIntVec Vector{JFInt}
typealias JFFltMat Matrix{JFFlt}
typealias JFIntMat Matrix{JFInt}
typealias JFMat{T<:Number} Matrix{T}
typealias JFVec{T<:Number} Vector{T}
export JFInt, JFFlt, JFCplxFlt, JFFltVec, JFIntVec, JFFltMat, JFIntMat, JFMat, JFVec

# Type for the model-data packaging system (used by all JFinEALE algorithms)
typealias ModelDataDictionary Dict{ASCIIString, Any}
export ModelDataDictionary

# Add key-value pairs to an existing model-data dictionary.
function dadd!(d::ModelDataDictionary; args...)
    for a in args
        key,val=a
        d[string(key)]=val;
    end
    return d
end
export dadd!

# Make model-data dictionary.
dmake(;args...)=dadd!(ModelDataDictionary(); args...)
export dmake

end

