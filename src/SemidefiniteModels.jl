module SemidefiniteModels

using MathProgBase
importall MathProgBase.SolverInterface

include("SD.jl")

include("sd_to_conic.jl")
include("conic_sdpa.jl")

end # module
