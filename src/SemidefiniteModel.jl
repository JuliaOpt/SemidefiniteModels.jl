module SemidefiniteModel

using MathProgBase
importall MathProgBase.SolverInterface

include("SD.jl")

include("sd_to_conic.jl")

end # module
