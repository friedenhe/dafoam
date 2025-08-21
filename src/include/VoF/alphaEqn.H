{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha() * mag(phi / mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        forAll(phic, faceI)
        {
            phic[faceI] *= (1.0 - icAlpha);
        }
        forAll(phic.boundaryField(), patchI)
        {
            forAll(phic.boundaryField()[patchI], faceI)
            {
                phic.boundaryFieldRef()[patchI][faceI] *= (1.0 - icAlpha);
            }
        }
        phic += (mixture.cAlpha() * icAlpha) * fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution
    if (scAlpha > 0)
    {
        phic +=
            scAlpha * mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }

    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    if (MULESCorr)
    {
#include "alphaSuSp.H"

        fvScalarMatrix alpha1Eqn(
            (
                fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1))
                + fv::gaussConvectionScheme<scalar>(
                      mesh,
                      phiCN,
                      upwind<scalar>(mesh, phiCN))
                      .fvmDiv(phiCN, alpha1)
            // - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
            //           + fvc::div(phiCN), alpha1)
            == Su + fvm::Sp(Sp + divU, alpha1));

        SolverPerformance<scalar> solverAlpha = alpha1Eqn.solve();
        DAUtility::primalResidualControl(solverAlpha, pimplePrintToScreen, "alpha", daGlobalVarPtr_->primalMaxRes);

        if (pimplePrintToScreen)
        {
            Info << alpha1.name() << " Mean: "
                 << alpha1.weightedAverage(mesh.Vsc()).value()
                 << "  Min: " << gMin(alpha1)
                 << "  Max: " << gMax(alpha1)
                 << endl;
        }

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    for (int aCorr = 0; aCorr < nAlphaCorr; aCorr++)
    {
#include "alphaSuSp.H"

        surfaceScalarField phir(phic * mixture.nHatf());

        tmp<surfaceScalarField> talphaPhi1Un(
            fvc::flux(
                phiCN(),
                alpha1,
                alphaScheme)
            + fvc::flux(
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme));

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField alpha10("alpha10", alpha1);

            MULESDF::correct(
                geometricOneField(),
                alpha1,
                talphaPhi1Un(),
                talphaPhi1Corr.ref(),
                Sp,
                (-Sp * alpha1)(),
                1,
                0);

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5 * alpha1 + 0.5 * alpha10;
                alphaPhi10 += 0.5 * talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un;

            MULESDF::explicitSolve(
                geometricOneField(),
                alpha1,
                phiCN,
                alphaPhi10,
                Sp,
                (Su + divU * min(alpha1(), scalar(1)))(),
                1,
                0);
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

#include "rhofs.H"

    rhoPhi = alphaPhi10 * (rho1f - rho2f) + phiCN * rho2f;

    if (pimplePrintToScreen)
    {
        Info << alpha1.name() << " Mean: "
             << alpha1.weightedAverage(mesh.Vsc()).value()
             << "  Min: " << gMin(alpha1)
             << "  Max: " << gMax(alpha1)
             << endl;
    }
}
