using SemidefiniteModels
using Base.Test

include("sdinterface.jl")
using CSDP
sdtest(CSDP.CSDPSolver(), duals=true)
