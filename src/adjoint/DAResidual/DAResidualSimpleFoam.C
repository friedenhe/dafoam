/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAResidualSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualSimpleFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualSimpleFoam::DAResidualSimpleFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
      setResidualClassMemberVector(U, dimensionSet(0, 1, -2, 0, 0, 0, 0)),
      setResidualClassMemberScalar(p, dimensionSet(0, 0, -1, 0, 0, 0, 0)),
      setResidualClassMemberPhi(phi),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create simpleControl
      simple_(const_cast<fvMesh&>(mesh))
{
}

void DAResidualSimpleFoam::calcResiduals(const dictionary& options)
{
    // We dont support MRF and fvOptions so all the related lines are commented
    // out for now

    // ******** U Residuals **********
    // copied and modified from UEqn.H

    word divUScheme = "div(phi,U)";

    label isPC = options.getLabel("isPC");

    if (isPC) 
    {
        divUScheme = "div(pc)";
    }

    tmp<fvVectorMatrix> tUEqn(
        fvm::div(phi_, U_, divUScheme)
        + daTurb_.divDevReff(U_));
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    URes_ = (UEqn & U_) + fvc::grad(p_);
    normalizeResiduals(URes);

    // ******** p Residuals **********
    // copied and modified from pEqn.H
    // NOTE manually set pRefCell and pRefValue
    label pRefCell = 0;
    scalar pRefValue = 0.0;

    volScalarField rAU(1.0 / UEqn.A());
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));
    //***************** NOTE *******************
    // we should not use the constrainHbyA function above since it
    // will degrade the accuracy of shape derivatives. Basically, we should
    // not constrain any variable because it will create discontinuity
    volVectorField HbyA("HbyA", U_);
    HbyA = rAU * UEqn.H();

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    adjustPhi(phiHbyA, U_, p_);

    tmp<volScalarField> rAtU(rAU);

    if (simple_.consistent())
    {
        rAtU = 1.0 / (1.0 / rAU - UEqn.H1());
        phiHbyA += fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p_) * mesh_.magSf();
        HbyA -= (rAU - rAtU()) * fvc::grad(p_);
    }

    tUEqn.clear();

    fvScalarMatrix pEqn(
        fvm::laplacian(rAtU(), p_)
        == fvc::div(phiHbyA));
    pEqn.setReference(pRefCell, pRefValue);

    pRes_ = pEqn & p_;
    normalizeResiduals(pRes);

    // ******** phi Residuals **********
    // copied and modified from pEqn.H
    phiRes_ = phiHbyA - pEqn.flux() - phi_;
    // need to normalize phiRes
    normalizePhiResiduals(phiRes);
}

void DAResidualSimpleFoam::updateIntermediateVariables()
{
}

/// update the boundary condition for all the states in the selected solver
void DAResidualSimpleFoam::correctBoundaryConditions()
{
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
