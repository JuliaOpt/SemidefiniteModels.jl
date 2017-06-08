# Methods for the Semidefinite interface

@compat abstract type AbstractSDModel <: AbstractMathProgModel end
export AbstractSDModel

MathProgBase.SolverInterface.@define_interface begin
    SDModel
    setconstrB!
    setconstrentry!
    setobjentry!
end
