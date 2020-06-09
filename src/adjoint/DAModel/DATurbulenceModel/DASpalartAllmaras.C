/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASpalartAllmaras.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASpalartAllmaras, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DASpalartAllmaras, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASpalartAllmaras::DASpalartAllmaras(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
      // SA parameters
      sigmaNut_(dimensioned<scalar>::lookupOrAddToDict(
          "sigmaNut",
          this->coeffDict_,
          0.66666)),
      kappa_(dimensioned<scalar>::lookupOrAddToDict(
          "kappa",
          this->coeffDict_,
          0.41)),

      Cb1_(dimensioned<scalar>::lookupOrAddToDict(
          "Cb1",
          this->coeffDict_,
          0.1355)),
      Cb2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cb2",
          this->coeffDict_,
          0.622)),
      Cw1_(Cb1_ / sqr(kappa_) + (1.0 + Cb2_) / sigmaNut_),
      Cw2_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw2",
          this->coeffDict_,
          0.3)),
      Cw3_(dimensioned<scalar>::lookupOrAddToDict(
          "Cw3",
          this->coeffDict_,
          2.0)),
      Cv1_(dimensioned<scalar>::lookupOrAddToDict(
          "Cv1",
          this->coeffDict_,
          7.1)),
      Cs_(dimensioned<scalar>::lookupOrAddToDict(
          "Cs",
          this->coeffDict_,
          0.3)),

      // Augmented variables
      nuTilda_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("nuTilda"))),
      nuTildaRes_(
          IOobject(
              "nuTildaRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
#ifdef CompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("nuTildaRes", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0),
#endif
          zeroGradientFvPatchScalarField::typeName),
      nuTildaResPartDeriv_(
          IOobject(
              "nuTildaResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          nuTildaRes_),
      y_(mesh.thisDb().lookupObject<volScalarField>("yWall"))
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// SA member functions
tmp<volScalarField> DASpalartAllmaras::chi() const
{
    return nuTilda_ / this->nu();
}

tmp<volScalarField> DASpalartAllmaras::fv1(
    const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3 / (chi3 + pow3(Cv1_));
}

tmp<volScalarField> DASpalartAllmaras::fv2(
    const volScalarField& chi,
    const volScalarField& fv1) const
{
    return 1.0 - chi / (1.0 + chi * fv1);
}

tmp<volScalarField> DASpalartAllmaras::Stilda(
    const volScalarField& chi,
    const volScalarField& fv1) const
{
    volScalarField Omega(::sqrt(2.0) * mag(skew(fvc::grad(U_))));

    return (
        max(
            Omega
                + fv2(chi, fv1) * nuTilda_ / sqr(kappa_ * y_),
            Cs_ * Omega));
}

tmp<volScalarField> DASpalartAllmaras::fw(
    const volScalarField& Stilda) const
{
    volScalarField r(
        min(
            nuTilda_
                / (max(
                       Stilda,
                       dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
                   * sqr(kappa_ * y_)),
            scalar(10.0)));
    r.boundaryFieldRef() == 0.0;

    const volScalarField g(r + Cw2_ * (pow6(r) - r));

    return g * pow((1.0 + pow6(Cw3_)) / (pow6(g) + pow6(Cw3_)), 1.0 / 6.0);
}

tmp<volScalarField> DASpalartAllmaras::DnuTildaEff() const
{
    return tmp<volScalarField>(
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu()) / sigmaNut_));
}

// Augmented functions
void DASpalartAllmaras::correctModelStates(wordList& modelStates) const
{
    // replace nut with nuTilda
    forAll(modelStates, idxI)
    {
        word stateName = modelStates[idxI];
        if (stateName == "nut")
        {
            modelStates[idxI] = "nuTilda";
        }
    }
}

/// update nut based on other turbulence variables and update the BCs
void DASpalartAllmaras::correctNut()
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    nut_ = nuTilda_ * fv1;

    nut_.correctBoundaryConditions();

    return;
}

/// update turbulence variable boundary values
void DASpalartAllmaras::correctBoundaryConditions()
{
    // correct the BCs for the perturbed fields
    nuTilda_.correctBoundaryConditions();
}

void DASpalartAllmaras::updateIntermediateVariables(const dictionary& options)
{
    // update nut based on nuTilda
    // Note: we need to update nut and its BC since we may have perturbed other turbulence vars
    // that affect the nut values
    this->correctNut();
}

void DASpalartAllmaras::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    // update the original variable connectivity for the adjoint state residuals in stateCon
    // For SA model just replace nut with nuTilda
    forAll(stateCon, idxI)
    {
        forAll(stateCon[idxI], idxJ)
        {
            word conStateName = stateCon[idxI][idxJ];
            if (conStateName == "nut")
            {
                stateCon[idxI][idxJ] = "nuTilda";
            }
        }
    }
}

void DASpalartAllmaras::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    // add the SA model residual connectivity to stateCon

    word pName;

    if (mesh_.thisDb().foundObject<volScalarField>("p"))
    {
        pName = "p";
    }
    else if (mesh_.thisDb().foundObject<volScalarField>("p_rgh"))
    {
        pName = "p_rgh";
    }
    else
    {
        FatalErrorIn(
            "Neither p nor p_rgh was found in mesh.thisDb()!"
            "addModelResidualCon failed to setup turbulence residuals!")
            << exit(FatalError);
    }

#ifdef IncompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "nuTilda", "phi"}, // lv0
            {"U", "nuTilda"}, // lv1
            {"nuTilda"} // lv2
        });
#endif

#ifdef CompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "T", pName, "nuTilda", "phi"}, // lv0
            {"U", "T", pName, "nuTilda"}, // lv1
            {"T", pName, "nuTilda"} // lv2
        });
#endif
}

/// solve the residual equations and update the state
void DASpalartAllmaras::correct()
{
    solveTurbState_ = 1;
    dictionary dummyOptions;
    this->calcResiduals(dummyOptions);

    solveTurbState_ = 0;
}

void DASpalartAllmaras::calcResiduals(const dictionary& options)
{
    // Copy and modify based on the "correct" function

    word divNuTildaScheme = "div(phi,nuTilda)";

    if (!solveTurbState_)
    {
        // we need to bound nuTilda before computing residuals
        // this will avoid having NaN residuals
        DAUtility::boundVar(allOptions_, nuTilda_);
    }

    //eddyViscosity<RASModelAugmented<BasicTurbulenceModel> >::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField Stilda(this->Stilda(chi, fv1));

    tmp<fvScalarMatrix> nuTildaEqn(
        fvm::ddt(phase_, rho_, nuTilda_)
            + fvm::div(phaseRhoPhi_, nuTilda_, divNuTildaScheme)
            - fvm::laplacian(phase_ * rho_ * DnuTildaEff(), nuTilda_)
            - Cb2_ / sigmaNut_ * phase_ * rho_ * magSqr(fvc::grad(nuTilda_))
        == Cb1_ * phase_ * rho_ * Stilda * nuTilda_
            - fvm::Sp(Cw1_ * phase_ * rho_ * fw(Stilda) * nuTilda_ / sqr(y_), nuTilda_));

    nuTildaEqn.ref().relax();

    if (solveTurbState_)
    {
        const scalar& deltaT = mesh_.time().deltaT().value();
        const scalar t = mesh_.time().timeOutputValue();
        label nSolverIters = round(t / deltaT);

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverNuTilda = solve(nuTildaEqn);
        if (nSolverIters % 100 == 0 || nSolverIters == 1)
        {
            Info << "nuTilda Initial residual: " << solverNuTilda.initialResidual() << endl
                 << "          Final residual: " << solverNuTilda.finalResidual() << endl;
        }

        DAUtility::boundVar(allOptions_, nuTilda_);
        nuTilda_.correctBoundaryConditions();

        // NOTE: in the original SA, it is correctNut(fv1) and fv1 is not
        // updated based on the latest nuTilda. We use correctNut which
        // recompute fv1 with the latest nuTilda
        this->correctNut();
    }
    else
    {
        // calculate residuals
        nuTildaRes_ = nuTildaEqn.ref() & nuTilda_;
    }

    return;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
