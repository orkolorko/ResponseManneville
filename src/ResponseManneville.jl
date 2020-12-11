module ResponseManneville

import IntervalArithmetic: Interval
export Interval

abstract type NormKind end
struct Linf <: NormKind end

include("DynamicDefinition.jl")
include("BasisDefinition.jl")
include("Contractors.jl")
include("GenericAssembler.jl")
include("C2Basis.jl")
include("InducedLSV.jl")


using .C2BasisDefinition, .InducedLSVMapDefinition
export ApproxInducedLSV, C2Basis

end