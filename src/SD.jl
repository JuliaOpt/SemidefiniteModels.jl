# Methods for the Semidefinite interface

abstract AbstractSDModel <: AbstractMathProgModel
export AbstractSDModel

MathProgBase.SolverInterface.@define_interface begin
    SDModel
    setconstrB!
    setconstrentry!
    setobjentry!
end
