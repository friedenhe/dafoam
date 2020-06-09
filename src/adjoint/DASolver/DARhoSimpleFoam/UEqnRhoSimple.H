// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn(
    fvm::div(phi, U)
    + turbulencePtr_->divDevRhoReff(U));
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

// get the solver performance info such as initial
// and final residuals
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));

if (nSolverIters % 100 == 0 || nSolverIters == 1)
{
    Info << "U Initial residual: " << solverU.initialResidual() << endl
         << "    Final residual: " << solverU.finalResidual() << endl;
}

// bound U
DAUtility::boundVar(allOptions, U);
U.correctBoundaryConditions();
