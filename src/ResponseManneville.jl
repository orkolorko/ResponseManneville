module ResponseManneville

import IntervalArithmetic: Interval
export Interval

abstract type NormKind end

struct L1 <: NormKind end
struct Linf <: NormKind end
struct Lipschitz <: NormKind end
struct TotalVariation <: NormKind end
struct ℓ1 <: NormKind end
struct ℓinf <: NormKind end

include("DynamicDefinition.jl")
include("BasisDefinition.jl")
include("Contractors.jl")
include("GenericAssembler.jl")
include("GenericEstimate.jl")
include("PwDynamicDefinition.jl")
include("Norms.jl")
include("NormsOfPowers.jl")
include("UlamBasis.jl")
include("C2Basis.jl")
include("InducedLSV.jl")
include("ContractionC1.jl")



using .C2BasisDefinition, .InducedLSVMapDefinition
export ApproxInducedLSV, C2Basis, C1, Ulam, assemble, DiscretizedOperator

end