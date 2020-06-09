volScalarField rAU(1.0 / UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho* rAU));
//volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
//***************** NOTE *******************
// we should not use the constrainHbyA function above since it
// will degrade the accuracy of shape derivatives. Basically, we should
// not constrain any variable because it will create discontinuity
volVectorField HbyA("HbyA", U);
HbyA = rAU * UEqn.H();
tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho) * fvc::flux(HbyA));

if (simple.transonic())
{
    surfaceScalarField phid(
        "phid",
        (fvc::interpolate(psi) / fvc::interpolate(rho)) * phiHbyA);

    phiHbyA -= fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(
            fvc::div(phiHbyA)
            + fvm::div(phid, p)
            - fvm::laplacian(rhorAUf, p));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.setReference(
            pressureControl.refCell(),
            pressureControl.refValue());

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverP = pEqn.solve();

        if (nSolverIters % 100 == 0 || nSolverIters == 1)
        {
            Info << "p Initial residual: " << solverP.initialResidual() << endl
                 << "    Final residual: " << solverP.finalResidual() << endl;
        }

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    closedVolume = adjustPhi(phiHbyA, U, p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(
            fvc::div(phiHbyA)
            - fvm::laplacian(rhorAUf, p));

        pEqn.setReference(
            pressureControl.refCell(),
            pressureControl.refValue());

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

if (nSolverIters % 100 == 0 || nSolverIters == 1)
{
#include "continuityErrsPython.H"
}

// Explicitly relax pressure for momentum corrector
p.relax();

// bound p
DAUtility::boundVar(allOptions, p);

U = HbyA - rAU * fvc::grad(p);
// bound U
DAUtility::boundVar(allOptions, U);
U.correctBoundaryConditions();

bool pLimited = pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass_ - fvc::domainIntegrate(psi * p))
        / fvc::domainIntegrate(psi);
}

if (pLimited || closedVolume)
{
    p.correctBoundaryConditions();
}

rho = thermo.rho();

if (!simple.transonic())
{
    rho.relax();
}

// bound rho
DAUtility::boundVar(allOptions, rho);
