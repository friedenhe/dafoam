/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAResidualRhoSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualRhoSimpleFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualRhoSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualRhoSimpleFoam::DAResidualRhoSimpleFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
      setResidualClassMemberVector(U, dimensionSet(1, -2, -2, 0, 0, 0, 0)),
      setResidualClassMemberScalar(p, dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      setResidualClassMemberScalar(T, dimensionSet(1, -1, -3, 0, 0, 0, 0)),
      setResidualClassMemberPhi(phi),
      thermo_(const_cast<fluidThermo&>(daModel.getThermo())),
      he_(thermo_.he()),
      rho_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("rho"))),
      alphat_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("alphat"))),
      psi_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("thermo:psi"))),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create simpleControl
      simple_(const_cast<fvMesh&>(mesh)),
      pressureControl_(p_, rho_, simple_.dict())
{

    // get molWeight and Cp from thermophysicalProperties
    const IOdictionary& thermoDict = mesh.thisDb().lookupObject<IOdictionary>("thermophysicalProperties");
    dictionary mixSubDict = thermoDict.subDict("mixture");
    dictionary specieSubDict = mixSubDict.subDict("specie");
    molWeight_ = specieSubDict.getScalar("molWeight");
    dictionary thermodynamicsSubDict = mixSubDict.subDict("thermodynamics");
    Cp_ = thermodynamicsSubDict.getScalar("Cp");

    // NOTE: for compressible flow, Prt is defined in RAS-turbulenceProperties
    // see EddyDiffusivity.C for reference
    const IOdictionary& turbDict = mesh.thisDb().lookupObject<IOdictionary>("turbulenceProperties");
    dictionary rasSubDict = turbDict.subDict("RAS");
    Prt_ = rasSubDict.getScalar("Prt");
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualRhoSimpleFoam::calcResiduals(const dictionary& options)
{
    // We dont support MRF and fvOptions so all the related lines are commented
    // out for now

    label isPC = options.getLabel("isPC");

    word divUScheme = "div(phi,U)";
    word divHEScheme = "div(phi,e)";
    word divPhidPScheme = "div(phid,p)";

    if (he_.name() == "h")
    {
        divHEScheme = "div(phi,h)";
    }

    if (isPC)
    {
        divUScheme = "div(pc)";
        divHEScheme = "div(pc)";
        divPhidPScheme = "div(pc)";
    }

    // ******** U Residuals **********
    // copied and modified from UEqn.H

    tmp<fvVectorMatrix> tUEqn(
        fvm::div(phi_, U_, divUScheme)
        + daTurb_.divDevRhoReff(U_));
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    URes_ = (UEqn & U_) + fvc::grad(p_);
    normalizeResiduals(URes);

    // ******** e Residuals **********
    // copied and modified from EEqn.H
    volScalarField alphaEff("alphaEff", thermo_.alphaEff(alphat_));

    fvScalarMatrix EEqn(
        fvm::div(phi_, he_, divHEScheme)
        + (he_.name() == "e"
               ? fvc::div(phi_, volScalarField("Ekp", 0.5 * magSqr(U_) + p_ / rho_))
               : fvc::div(phi_, volScalarField("K", 0.5 * magSqr(U_))))
        - fvm::laplacian(alphaEff, he_));

    EEqn.relax();

    TRes_ = EEqn & he_;
    normalizeResiduals(TRes);

    // ******** p and phi Residuals **********
    // copied and modified from pEqn.H
    volScalarField rAU(1.0 / UEqn.A());
    volScalarField rAtU(1.0 / (1.0 / rAU - UEqn.H1()));
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    //***************** NOTE *******************
    // we should not use the constrainHbyA function above since it
    // will degrade the accuracy of shape derivatives. Basically, we should
    // not constrain any variable because it will create discontinuity
    volVectorField HbyA("HbyA", U_);
    HbyA = rAU * UEqn.H();
    tUEqn.clear();

    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho_) * fvc::flux(HbyA));

    volScalarField rhorAtU("rhorAtU", rho_ * rAtU);

    // Update the pressure BCs to ensure flux consistency
    // constrainPressure(p_, rho_, U_, phiHbyA, rhorAtU, this->MRF_);

    // NOTE: we don't support transonic = true

    adjustPhi(phiHbyA, U_, p_);

    phiHbyA += fvc::interpolate(rho_ * (rAtU - rAU)) * fvc::snGrad(p_) * mesh_.magSf();
    HbyA -= (rAU - rAtU) * fvc::grad(p_);

    fvScalarMatrix pEqn(
        fvc::div(phiHbyA)
        - fvm::laplacian(rhorAtU, p_));

    pEqn.setReference(pressureControl_.refCell(), pressureControl_.refValue());

    pRes_ = pEqn & p_;
    normalizeResiduals(pRes);

    // ******** phi Residuals **********
    // copied and modified from pEqn.H
    phiRes_ = phiHbyA + pEqn.flux() - phi_;
    normalizePhiResiduals(phiRes);
}

void DAResidualRhoSimpleFoam::updateIntermediateVariables()
{
    // ********************** NOTE *****************
    // we assume hePsiThermo
    // TODO: need to do this using built-in openfoam functions.

    // we need to:
    // 1, update psi based on T, psi=1/(R*T)
    // 2, update rho based on p and psi, rho=psi*p
    // 3, update E based on T, p and rho, E=Cp*T-p/rho
    // 4, update alphat
    // 5, update velocity boundary based on MRF

    scalar RR = Foam::constant::thermodynamic::RR; // 8314.4700665  gas constant in OpenFOAM
    dimensionedScalar R(
        "R",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        RR / molWeight_);
    psi_ = 1 / T_ / R;
    rho_ = psi_ * p_;

    forAll(he_, idxI)
    {
        if (he_.name() == "e")
        {
            he_[idxI] = Cp_ * T_[idxI] - p_[idxI] / rho_[idxI];
        }
        else
        {
            he_[idxI] = Cp_ * T_[idxI];
        }
    }
    he_.correctBoundaryConditions();

    // NOTE: for compressible flow, Prt is defined in RAS-turbulenceProperties
    // see EddyDiffusivity.C for reference
    dimensionedScalar Prt1(
        "Prt1",
        dimless,
        Prt_);

    alphat_ = rho_ * daTurb_.getNut() / Prt1;
    alphat_.correctBoundaryConditions();

}

/// update the boundary condition for all the states in the selected solver
void DAResidualRhoSimpleFoam::correctBoundaryConditions()
{
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
}

} // End namespace Foam

// ************************************************************************* //
