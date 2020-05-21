// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn(
    fvm::div(phi, U)
    + turbulencePtr_->divDevRhoReff(U));
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

solve(UEqn == -fvc::grad(p));

// bound U
daUtilPtr_->boundVar(allOptions, U);
U.correctBoundaryConditions();
