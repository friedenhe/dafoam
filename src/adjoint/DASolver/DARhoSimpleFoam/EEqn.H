{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn(
        fvm::div(phi, he)
        + (he.name() == "e"
               ? fvc::div(phi, volScalarField("Ekp", 0.5 * magSqr(U) + p / rho))
               : fvc::div(phi, volScalarField("K", 0.5 * magSqr(U))))
        - fvm::laplacian(turbulencePtr_->alphaEff(), he));

    EEqn.relax();

    EEqn.solve();

    // bound he
    daUtilPtr_->boundVar(allOptions, he);

    thermo.correct();
}
