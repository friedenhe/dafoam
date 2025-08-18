/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAResidualInterIsoFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualInterIsoFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualInterIsoFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualInterIsoFoam::DAResidualInterIsoFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
      setResidualClassMemberVector(U, dimensionSet(1, -2, -2, 0, 0, 0, 0)),
      setResidualClassMemberScalar(p_rgh, dimensionSet(0, 0, -1, 0, 0, 0, 0)),
      alpha1_(mesh_.thisDb().lookupObjectRef<volScalarField>("alpha.water")),
      alpha1Res_(
          IOobject(
              "alpha.waterRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimensionedScalar("alpha.waterRes", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchField<scalar>::typeName),
      setResidualClassMemberPhi(phi),
      rho_(mesh_.thisDb().lookupObjectRef<volScalarField>("rho")),
      rhoPhi_(mesh_.thisDb().lookupObjectRef<surfaceScalarField>("rhoPhi")),
      gh_(mesh_.thisDb().lookupObjectRef<volScalarField>("gh")),
      ghf_(mesh_.thisDb().lookupObjectRef<surfaceScalarField>("ghf")),
      advector_(mesh_.thisDb().lookupObjectRef<isoAdvectionDF>("isoAdvectionDF")),
      mixture_(mesh_.thisDb().lookupObjectRef<immiscibleIncompressibleTwoPhaseMixture>("transportProperties")),
      turbulence_(mesh_.thisDb().lookupObjectRef<incompressible::turbulenceModel>("turbulenceProperties")),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create simpleControl
      pimple_(const_cast<fvMesh&>(mesh))
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualInterIsoFoam::clear()
{
    /*
    Description:
        Clear all members to avoid memory leak because we will initalize 
        multiple objects of DAResidual. Here we need to delete all members
        in the parent and child classes
    */

    URes_.clear();
    p_rghRes_.clear();
    phiRes_.clear();
    alpha1Res_.clear();
}

void DAResidualInterIsoFoam::calcResiduals(const dictionary& options)
{
    /*
    Description:
        This is the function to compute residuals.
    
    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,U)

        p_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        URes_, pRes_, phiRes_: residual field variables
    */

    // ******** alpha1 Residuals **********
    // copied and modified from alphaEqn.H

    // Initialising dVf with upwind values
    // i.e. phi[facei]*alpha1[upwindCell[facei]]*dt
    surfaceScalarField dVf = upwind<scalar>(mesh_, phi_).flux(alpha1_) * mesh_.time().deltaT();

    advector_.timeIntegratedFlux(dVf);

    //advector_.limitFluxes();
    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();
    rhoPhi_ = (rho1 - rho2) * dVf / mesh_.time().deltaT() + rho2 * phi_;

    alpha1Res_ = alpha1_.oldTime() - fvc::surfaceIntegrate(dVf) - alpha1_;
    normalizeResiduals(alpha1Res);

    // ******** U Residuals **********
    // copied and modified from UEqn.H

    word divUScheme = "div(rhoPhi,U)";

    label isPC = options.getLabel("isPC");

    if (isPC)
    {
        divUScheme = "div(pc)";
    }

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
        + fvm::div(rhoPhi_, U_, divUScheme)
        + turbulence_.divDevRhoReff(rho_, U_));

    // NOTE: we need to call UEqn.relax here because it does some BC treatment, but we need to
    // force the relaxation factor to be 1.0 because the last pimple loop does not use relaxation
    UEqn.relax(1.0);

    URes_ = (UEqn & U_)
        - fvc::reconstruct((mixture_.surfaceTensionForce()
                            - ghf_ * fvc::snGrad(rho_) - fvc::snGrad(p_rgh_))
                           * mesh_.magSf());
    normalizeResiduals(URes);

    // ******** p Residuals **********
    // copied and modified from pEqn.H
    // NOTE manually set pRefCell and pRefValue
    label pRefCell = 0;
    scalar pRefValue = 0.0;

    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    //***************** NOTE *******************
    // constrainHbyA has been used since OpenFOAM-v1606; however, it may degrade the accuracy of derivatives
    // because constraining variables will create discontinuity. Here we have a option to use the old
    // implementation in OpenFOAM-3.0+ and before (no constraint for HbyA)
    autoPtr<volVectorField> HbyAPtr = nullptr;
    label useConstrainHbyA = daOption_.getOption<label>("useConstrainHbyA");
    if (useConstrainHbyA)
    {
        HbyAPtr.reset(new volVectorField(constrainHbyA(rAU * UEqn.H(), U_, p_rgh_)));
    }
    else
    {
        HbyAPtr.reset(new volVectorField("HbyA", U_));
        HbyAPtr() = rAU * UEqn.H();
    }
    volVectorField& HbyA = HbyAPtr();

    surfaceScalarField phiHbyA(
        "phiHbyA",
        fvc::flux(HbyA));

    surfaceScalarField phig(
        (
            mixture_.surfaceTensionForce()
            - ghf_ * fvc::snGrad(rho_))
        * rAUf * mesh_.magSf());

    phiHbyA += phig;

    fvScalarMatrix p_rghEqn(
        fvm::laplacian(rAUf, p_rgh_)
        == fvc::div(phiHbyA));

    p_rghEqn.setReference(pRefCell, pRefValue);

    p_rghRes_ = p_rghEqn & p_rgh_;
    normalizeResiduals(p_rghRes);

    // ******** phi Residuals **********
    // copied and modified from pEqn.H
    phiRes_ = phiHbyA - p_rghEqn.flux() - phi_;
    // need to normalize phiRes
    normalizePhiResiduals(phiRes);
}

void DAResidualInterIsoFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /* 
    Description:
        Calculate the diagonal block of the preconditioner matrix dRdWTPC using the fvMatrix
    */
}

void DAResidualInterIsoFoam::updateIntermediateVariables()
{
    /* 
    Description:
        Update the intermediate variables that depend on the state variables
    */

    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();

    volScalarField& alpha1(mixture_.alpha1());
    volScalarField& alpha2(mixture_.alpha2());
    alpha2 = 1.0 - alpha1;
    mixture_.correct();
    rho_ == alpha1* rho1 + alpha2* rho2;

    volScalarField& p = mesh_.thisDb().lookupObjectRef<volScalarField>("p");
    p == p_rgh_ + rho_* gh_;
}

void DAResidualInterIsoFoam::correctBoundaryConditions()
{
    /* 
    Description:
        Update the boundary condition for all the states in the selected solver
    */

    U_.correctBoundaryConditions();
    p_rgh_.correctBoundaryConditions();
    alpha1_.correctBoundaryConditions();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
