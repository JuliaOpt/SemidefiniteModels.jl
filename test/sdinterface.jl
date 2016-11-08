using Base.Test
using SemidefiniteModel
using MathProgBase.SolverInterface

function sdtest(solver::MathProgBase.AbstractMathProgSolver; duals=false, tol=1e-6)
    m = ConicModel(solver)
    loadproblem!(m, "prob.dat-s")
    optimize!(m)
    @test status(m) == :Optimal
    @test isapprox(getobjval(m), 2.75)
    @test norm(getsolution(m) - [.75, 1.]) < tol
    if duals
        @test norm(getdual(m) - [0, 0, .125, .125*sqrt(2), .125, 2/3, 0, 0, 0, 0, 0]) < tol
        @test norm(getvardual(m)) < tol
    end
end
