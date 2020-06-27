{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn(
        fvm::div(phi, he)
        + (he.name() == "e"
               ? fvc::div(phi, volScalarField("Ekp", 0.5 * magSqr(U) + p / rho))
               : fvc::div(phi, volScalarField("K", 0.5 * magSqr(U))))
        - fvm::laplacian(turbulencePtr_->alphaEff(), he));

    EEqn.relax();

    // get the solver performance info such as initial
    // and final residuals
    SolverPerformance<scalar> solverE = EEqn.solve();

    if (nSolverIters % 100 == 0 || nSolverIters == 1)
    {
        Info << "he Initial residual: " << solverE.initialResidual() << endl
             << "     Final residual: " << solverE.finalResidual() << endl;
    }

    // bound he
    DAUtility::boundVar(allOptions, he);

    thermo.correct();
}
