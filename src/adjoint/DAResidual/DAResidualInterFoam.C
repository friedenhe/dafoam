/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAResidualInterFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualInterFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualInterFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualInterFoam::DAResidualInterFoam(
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
          dimensionedScalar("alpha.waterRes", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
          zeroGradientFvPatchField<scalar>::typeName),
      setResidualClassMemberPhi(phi),
      rho_(mesh_.thisDb().lookupObjectRef<volScalarField>("rho")),
      rhoPhi_(mesh_.thisDb().lookupObjectRef<surfaceScalarField>("rhoPhi")),
      gh_(mesh_.thisDb().lookupObjectRef<volScalarField>("gh")),
      ghf_(mesh_.thisDb().lookupObjectRef<surfaceScalarField>("ghf")),
      alphaPhi10_(mesh_.thisDb().lookupObjectRef<surfaceScalarField>("alphaPhi10")),
      mixture_(mesh_.thisDb().lookupObjectRef<immiscibleIncompressibleTwoPhaseMixture>("transportProperties")),
      turbulence_(mesh_.thisDb().lookupObjectRef<incompressible::turbulenceModel>("turbulenceProperties")),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create simpleControl
      pimple_(const_cast<fvMesh&>(mesh)),
      MULESDF_(mesh_.thisDb().lookupObjectRef<MULESDF>("MULESDF"))
{
    const dictionary& alphaControls = mesh.solverDict(alpha1_.name());

    MULESCorr_ = alphaControls.lookupOrDefault("MULESCorr", false);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualInterFoam::clear()
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

void DAResidualInterFoam::calcResiduals(const dictionary& options)
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

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    volScalarField& alpha1(mixture_.alpha1());
    volScalarField& alpha2(mixture_.alpha2());
    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();

    surfaceScalarField phic(mixture_.cAlpha() * mag(phi_ / mesh_.magSf()));

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

    tmp<surfaceScalarField> phiCN(phi_);

    if (MULESCorr_)
    {
#include "alphaSuSp.H"
        fvScalarMatrix alpha1Eqn(
            (
                fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha1))
                + fv::gaussConvectionScheme<scalar>(
                      mesh_,
                      phiCN,
                      upwind<scalar>(mesh_, phiCN))
                      .fvmDiv(phiCN, alpha1)
            // - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
            //           + fvc::div(phiCN), alpha1)
            == Su + fvm::Sp(Sp + divU, alpha1));
        alpha1Res_ = alpha1Eqn & alpha1_;
        normalizeResiduals(alpha1Res);
    }
    else
    {

        surfaceScalarField phir(phic * mixture_.nHatf());

        tmp<surfaceScalarField> talphaPhi1Un(
            fvc::flux(
                phiCN(),
                alpha1,
                alphaScheme)
            + fvc::flux(
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme));

        alphaPhi10_ = talphaPhi1Un;

        const scalar& rDeltaT = 1.0 / mesh_.time().deltaTValue();

        scalarField psiIf = alpha1;
        const scalarField& psi0 = alpha1.oldTime();

        forAll(psiIf, cellI)
        {
            psiIf[cellI] = 0.0;
        }
        fvc::surfaceIntegrate(psiIf, alphaPhi10_);

        forAll(alpha1Res_, cellI)
        {
            alpha1Res_[cellI] = ((rho_.oldTime()[cellI] * psi0[cellI] * rDeltaT - psiIf[cellI])
                                     / (rho_[cellI] * rDeltaT)
                                 - alpha1_[cellI])
                * rDeltaT / mesh_.V()[cellI];
        }
        normalizeResiduals(alpha1Res);
        // update rhoPhi here
        rhoPhi_ = alphaPhi10_ * (rho1 - rho2) + phiCN * rho2;
    }

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

    //adjustPhi(phiHbyA, U_, p_rgh_);

    surfaceScalarField phig(
        (
            mixture_.surfaceTensionForce()
            - ghf_ * fvc::snGrad(rho_))
        * rAUf * mesh_.magSf());

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    //constrainPressure(p_rgh_, U_, phiHbyA, rAUf);

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

void DAResidualInterFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /* 
    Description:
        Calculate the diagonal block of the preconditioner matrix dRdWTPC using the fvMatrix
    */

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    PetscScalar val;

    dictionary normStateDict = daOption_.getAllOptions().subDict("normalizeStates");
    wordList normResDict = daOption_.getOption<wordList>("normalizeResiduals");

    scalar UScaling = 1.0;
    if (normStateDict.found("U"))
    {
        UScaling = normStateDict.getScalar("U");
    }
    scalar UResScaling = 1.0;

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
        + fvm::div(rhoPhi_, U_, "div(pc)")
        + turbulence_.divDevRhoReff(rho_, U_));

    UEqn.relax(1.0);

    // set diag
    forAll(U_, cellI)
    {
        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[cellI];
        }
        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", cellI, i);
            PetscInt colI = rowI;
            scalarField D = UEqn.D();
            scalar val1 = D[cellI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
        }
    }

    // set lower/owner
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[neighbourCellI];
        }

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            scalar val1 = UEqn.lower()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    // set upper/neighbour
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[ownerCellI];
        }

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            scalar val1 = UEqn.upper()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
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

    // ********* p
    scalar pScaling = 1.0;
    if (normStateDict.found("p"))
    {
        pScaling = normStateDict.getScalar("p");
    }
    scalar pResScaling = 1.0;
    // set diag
    forAll(p_rgh_, cellI)
    {
        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[cellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p_rgh", cellI);
        PetscInt colI = rowI;
        scalarField D = p_rghEqn.D();
        scalar val1 = D[cellI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // set lower/owner
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[neighbourCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p_rgh", neighbourCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p_rgh", ownerCellI);
        scalar val1 = p_rghEqn.lower()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    // set upper/neighbour
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[ownerCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p_rgh", ownerCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p_rgh", neighbourCellI);
        scalar val1 = p_rghEqn.upper()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }
}

void DAResidualInterFoam::updateIntermediateVariables()
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

    // TODO this is a simplification for now.
    if (MULESCorr_)
    {
        rhoPhi_ = fvc::interpolate(rho_) * phi_;
    }

    volScalarField& p = mesh_.thisDb().lookupObjectRef<volScalarField>("p");
    p == p_rgh_ + rho_* gh_;
}

void DAResidualInterFoam::correctBoundaryConditions()
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
